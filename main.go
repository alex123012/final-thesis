package main

import (
	"bytes"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/alex123012/final-thesis/pkg/shell"
	"github.com/biogo/biogo/seq/linear"
)

const (
	bowtieIndexFolder        = "bowtie_indexes"
	adapterTrimmedFilePrefix = "trimmed_"
	resultsFolder            = "results"
)

var (
	GenomeVersion           = "GCF_000001405.40"
	GenomeAnnotationGffFile = "mirna/hsa.gff3"
	MirnaHairpinSeqsFile    = "mirna/hairpin.fa"
	MirnaMatureSeqsFile     = "mirna/mature.fa"
	ChrFileMapper           = "full_mapper.tsv"
	GenomeFile              string
	BowtieIndexPrefix       string
	Accession               string
	AdapterSeq              string

	fastqFile                   string
	accessionResultsFolder      string
	bowtieResultFile            string
	bamResultFile               string
	bamSortedResultFile         string
	vcfSortedFilteredResultFile string
	tsvSortedFilteredResultFile string
	resultTable                 string
	resultHeatmapImage          string
	resultBarplotImage          string

	reheaderedGenomeFile    string
	genomeAnnotationBedFile string
)

func parseArgs() error {

	flag.StringVar(&GenomeVersion, "genome-version", GenomeVersion, "genome version to download/extract files from.\nWon't be used if --genome-file exists.\noptional")
	flag.StringVar(&ChrFileMapper, "chr-file-mapper", ChrFileMapper, "tsv file with two columns: chr id and chr* name")
	flag.StringVar(&GenomeAnnotationGffFile, "annotation-file", GenomeAnnotationGffFile, "gff file containing genome regions annotations")
	flag.StringVar(&MirnaHairpinSeqsFile, "mirna-hairpin-file", MirnaHairpinSeqsFile, "file path to mirna hairpin seqs")
	flag.StringVar(&MirnaMatureSeqsFile, "mirna-mature-file", MirnaMatureSeqsFile, "file path to mature mirna seqs")

	flag.StringVar(&GenomeFile, "genome-file", GenomeFile, "genome file name, if not exists - will be extracted from downloaded NCBI archive\noptional")
	flag.StringVar(&BowtieIndexPrefix, "bowtie-index-prefix", BowtieIndexPrefix, "bowtie index prefix - with this prefix bowtie2-build will create\nindex files if they not exists.\noptional")
	flag.StringVar(&Accession, "accession", Accession, "Accession number of SRA experiment and folder name in `results` dir to save results in")
	flag.StringVar(&AdapterSeq, "adapter-seq", AdapterSeq, "if provided, cutadapt will be used to trim provided adapters seq")
	flag.Parse()

	if GenomeVersion == "" && GenomeFile == "" {
		return errors.New("provide genome version or genome file")
	}

	if GenomeFile == "" {
		GenomeFile = filepath.Join("genomes", GenomeVersion+".fasta")
	}

	if BowtieIndexPrefix == "" {
		BowtieIndexPrefix = GenomeVersion
	}

	if BowtieIndexPrefix == "" {
		return errors.New("provide bowtie index prefix")
	}

	if Accession == "" {
		return errors.New("provide SRA accession name")
	}
	if ChrFileMapper == "" {
		return errors.New("provide chr mapper file")
	}
	if GenomeAnnotationGffFile == "" {
		return errors.New("provide genome annotation gff file")
	}
	if MirnaHairpinSeqsFile == "" {
		return errors.New("provide mirna hairpin seqs file")
	}
	if MirnaMatureSeqsFile == "" {
		return errors.New("provide mirna mature seqs file")
	}

	BowtieIndexPrefix = filepath.Join(bowtieIndexFolder, BowtieIndexPrefix)
	fastqFile = filepath.Join("input", Accession+".fastq.gz")
	accessionResultsFolder = filepath.Join(resultsFolder, "samples_data", Accession)

	bowtieResultFile = filepath.Join(accessionResultsFolder, "result.sam")
	bamResultFile = filepath.Join(accessionResultsFolder, "result.bam")
	bamSortedResultFile = filepath.Join(accessionResultsFolder, "result.sorted.bam")
	vcfSortedFilteredResultFile = filepath.Join(accessionResultsFolder, "result.sorted.filtered.vcf")
	tsvSortedFilteredResultFile = filepath.Join(accessionResultsFolder, "result.sorted.filtered.annotated.tsv")
	resultTable = filepath.Join(resultsFolder, "tables/result-"+Accession+".tsv")
	resultHeatmapImage = filepath.Join(resultsFolder, "heatmaps/result-heatmap-"+Accession+".png")
	resultBarplotImage = filepath.Join(resultsFolder, "barplots/result-barplot-"+Accession+".png")

	reheaderedGenomeFile = removeExtFromFile(GenomeFile) + ".reheadered" + filepath.Ext(GenomeFile)

	genomeAnnotationBedFile = removeExtFromFile(GenomeAnnotationGffFile) + ".bed"
	return nil
}

func removeExtFromFile(f string) string {
	return strings.TrimSuffix(f, filepath.Ext(f))
}

func main() {
	// brew install samtools bcftools bedtools bedops
	if err := parseArgs(); err != nil {
		log.Fatalln(err)
	}

	log.Println("Running fasterq-dump")
	if err := runFasterqDump(Accession, fastqFile); err != nil {
		log.Fatalln(err)
	}

	if AdapterSeq != "" {
		log.Println("Running cutadapt")
		trimmedFastqFiles, err := runCutadapt(AdapterSeq, adapterTrimmedFilePrefix, accessionResultsFolder, fastqFile)
		if err != nil {
			log.Fatalln(err)
		}
		fastqFile = trimmedFastqFiles[0]
	}

	chrMapper, err := chrMapper(ChrFileMapper)
	if err != nil {
		log.Fatalln(err)
	}

	runConfigs := []struct {
		resultFile string
		condition  func() bool
		f          func() error
		logMessage string
	}{
		// Downloading step
		{
			resultFile: GenomeFile,
			f:          getGenomeByVersion(GenomeVersion, GenomeFile),
			logMessage: fmt.Sprintf("Downloading genome with version %s", GenomeVersion),
		},
		{
			resultFile: GenomeAnnotationGffFile,
			f:          getMirnaAnnotations(GenomeAnnotationGffFile),
			logMessage: fmt.Sprintf("Downloading miRNA annotations file %s", GenomeAnnotationGffFile),
		},
		{
			resultFile: MirnaMatureSeqsFile,
			f:          getMirnaSeqs("mature", MirnaMatureSeqsFile),
			logMessage: fmt.Sprintf("Downloading miRNA mature seqs %s", MirnaMatureSeqsFile),
		},
		{
			resultFile: MirnaHairpinSeqsFile,
			f:          getMirnaSeqs("hairpin", MirnaHairpinSeqsFile),
			logMessage: fmt.Sprintf("Downloading miRNA hairpin seqs %s", MirnaHairpinSeqsFile),
		},
		// Preparation step
		{
			resultFile: reheaderedGenomeFile,
			f:          reheaderGenomeFile(GenomeFile, reheaderedGenomeFile, chrMapper),
			logMessage: "Running genome file reheadering",
		},
		{
			resultFile: genomeAnnotationBedFile,
			f:          gffToBed(GenomeAnnotationGffFile, genomeAnnotationBedFile),
			logMessage: "Running gff2bed",
		},
		// Results step
		{
			resultFile: BowtieIndexPrefix + ".1.ebwt",
			f:          runBowtieBuild(reheaderedGenomeFile, BowtieIndexPrefix),
			logMessage: "Running bowtie-build",
		},
		{
			resultFile: bowtieResultFile,
			condition:  func() bool { return fileExists(bowtieResultFile + ".gz") },
			f:          runBowtie(BowtieIndexPrefix, bowtieResultFile, fastqFile),
			logMessage: "Running bowtie",
		},
		{
			resultFile: bamResultFile,
			f:          runSamtoolsViewBam(bowtieResultFile, bamResultFile),
			logMessage: "Running samtools view -b",
		},
		{
			resultFile: bowtieResultFile + ".gz",
			f:          gzipFile(bowtieResultFile),
			logMessage: "Saving bowtie result (SAM) file to gzip archive",
		},
		{
			resultFile: bamSortedResultFile,
			f:          runSamtoolsViewSort(bamResultFile, bamSortedResultFile),
			logMessage: "Running samtools sort",
		},
		{
			resultFile: bamSortedResultFile + ".bai",
			f:          runSamtoolsIndex(bamSortedResultFile),
			logMessage: "Running samtools index",
		},
		{
			resultFile: vcfSortedFilteredResultFile,
			f:          runBcftoolsMpileup(reheaderedGenomeFile, genomeAnnotationBedFile, bamSortedResultFile, vcfSortedFilteredResultFile),
			logMessage: "Running bcftools mpileup",
		},
		{
			resultFile: tsvSortedFilteredResultFile,
			f:          runBedtoolsIntersectMpileupAndGff(GenomeAnnotationGffFile, vcfSortedFilteredResultFile, tsvSortedFilteredResultFile),
			logMessage: "Running bedtools intersect",
		},
		{
			resultFile: resultTable,
			f:          generateResultTable(MirnaMatureSeqsFile, MirnaHairpinSeqsFile, reheaderedGenomeFile, tsvSortedFilteredResultFile, resultTable),
			logMessage: "Generating result tsv file for plotting",
		},
		{
			resultFile: resultHeatmapImage,
			f:          generateHeatmap(resultTable, resultHeatmapImage),
			logMessage: "Generating Heatmap image",
		},
		{
			resultFile: resultBarplotImage,
			f:          generateBarplot(resultTable, resultBarplotImage),
			logMessage: "Generating Barplot image",
		},
	}

	for _, cfg := range runConfigs {
		if fileExists(cfg.resultFile) || (cfg.condition != nil && cfg.condition()) {
			continue
		}

		log.Println(cfg.logMessage)
		if err := cfg.f(); err != nil {
			log.Fatalln(err)
		}
	}
	log.Println("Work is done")
}

func runFasterqDump(accession, resultFile string) error {
	if fileExists(resultFile) {
		return nil
	}

	if err := pathForFile(resultFile); err != nil {
		return err
	}

	fastqFileUngzipped := removeExtFromFile(resultFile)

	sr := shell.NewShellRunner(nil, nil, nil)
	if err := sr.RunShell("fasterq-dump", "--outfile", fastqFileUngzipped, "--mem", "1GB", "-x", accession); err != nil {
		return err
	}

	return gzipFile(fastqFileUngzipped)()
}

func runCutadapt(adapterSeq, resultPrefix, resultsFolder string, fastqFiles ...string) ([]string, error) {
	if err := mkdir(resultsFolder); err != nil {
		return nil, err
	}

	sr := shell.NewShellRunner(nil, nil, nil)
	trimmedFastqFiles := make([]string, 0, len(fastqFiles))
	for _, fqFile := range fastqFiles {
		resultFqFile := filepath.Join(resultsFolder, resultPrefix+filepath.Base(fqFile))
		trimmedFastqFiles = append(trimmedFastqFiles, resultFqFile)
		if fileExists(resultFqFile) {
			continue
		}
		if err := sr.RunShell("cutadapt", "-q", "15", "-e", "0.12", "-a", adapterSeq, "-m", "16", "--discard-untrimmed", fqFile, "--output", resultFqFile); err != nil {
			return nil, err
		}
	}
	return trimmedFastqFiles, nil
}

func getGenomeByVersion(genomeVersion, genomeFile string) func() error {
	return func() error {
		if genomeVersion == "" {
			return errors.New("genome file doesn't exists, please provide genome version")
		}

		if err := pathForFile(genomeFile); err != nil {
			return err
		}
		zipFile, err := downloadGenomeSequence(genomeVersion, genomeFile)
		if err != nil {
			return err
		}

		fileInArchive := fmt.Sprintf("ncbi_dataset/data/%s/%s_GRCh38.p14_genomic.fna", genomeVersion, genomeVersion)
		if err := extractFileFromZip(zipFile, fileInArchive, genomeFile); err != nil {
			return err
		}
		return nil
	}
}

func getMirnaAnnotations(annotationsFile string) func() error {
	return func() error {
		if err := pathForFile(annotationsFile); err != nil {
			return err
		}
		return downloadMirnaAnnotations(annotationsFile)
	}
}

func getMirnaSeqs(seqsType, filePath string) func() error {
	return func() error {
		if err := pathForFile(filePath); err != nil {
			return err
		}

		return downloadMirnaSeqs(seqsType, filePath)
	}
}

func reheaderGenomeFile(genomeFile, resultFile string, mapper map[string]string) func() error {
	return func() error {
		if err := pathForFile(resultFile); err != nil {
			return err
		}

		rf, err := os.Open(genomeFile)
		if err != nil {
			return err
		}
		defer rf.Close()

		wf, err := os.Create(resultFile)
		if err != nil {
			return err
		}
		defer wf.Close()

		sc := newDNAFastaReader(rf)
		for sc.Next() {
			s := sc.Seq().(*linear.Seq)
			s.SetName(mapper[s.ID])
			if _, err := wf.Write(seqToBytes(s)); err != nil {
				return err
			}
		}
		if err := sc.Error(); err != nil {
			return err
		}
		return nil
	}
}

func gffToBed(annotationFileGff, annotationFileBed string) func() error {
	return func() error {
		if err := pathForFile(annotationFileBed); err != nil {
			return err
		}

		rf, err := os.Open(annotationFileGff)
		if err != nil {
			return err
		}
		defer rf.Close()

		wf, err := os.Create(annotationFileBed)
		if err != nil {
			return err
		}
		defer wf.Close()

		sr := shell.NewShellRunner(rf, wf, nil)
		return sr.RunShell("gff2bed")
	}
}

func runBowtieBuild(genomeFile, bowtieIndexPrefix string) func() error {
	return func() error {
		if err := pathForFile(bowtieIndexPrefix); err != nil {
			return err
		}
		sr := shell.NewShellRunner(nil, nil, nil)
		return sr.RunShell("bowtie-build", genomeFile, bowtieIndexPrefix)
	}
}

func runBowtie(bowtieIndexPrefix, resultsFilename string, fastqFiles ...string) func() error {
	return func() error {
		if err := pathForFile(resultsFilename); err != nil {
			return err
		}

		sr := shell.NewShellRunner(nil, nil, nil)
		return sr.RunShell("bowtie", "-S", "-v", "1", "-a", "--best", "--strata", "-x", bowtieIndexPrefix, strings.Join(fastqFiles, ","), resultsFilename)
	}
}

func runSamtools(args ...string) error {
	sr := shell.NewShellRunner(nil, nil, nil)
	return sr.RunShell("samtools", args...)
}

func runSamtoolsViewBam(samFile, resultBamFile string) func() error {
	return func() error {
		if err := pathForFile(resultBamFile); err != nil {
			return err
		}
		if !fileExists(samFile) {
			sr := shell.NewShellRunner(nil, nil, nil)
			if err := sr.RunShell("gunzip", samFile+".gz"); err != nil {
				return err
			}
		}
		return runSamtools("view", "-bo", resultBamFile, samFile)
	}
}

func runSamtoolsViewSort(bamFile, resultSortedFile string) func() error {
	return func() error {
		if err := pathForFile(resultSortedFile); err != nil {
			return err
		}
		return runSamtools("sort", "-o", resultSortedFile, bamFile)
	}
}

func runSamtoolsIndex(bamFile string) func() error {
	return func() error {
		if err := pathForFile(bamFile); err != nil {
			return err
		}
		return runSamtools("index", bamFile)
	}
}

func runBcftoolsMpileup(genomeFile, annotationFile, bamFile, resultFile string) func() error {
	return func() error {
		if err := pathForFile(resultFile); err != nil {
			return err
		}

		outputfile, err := os.Create(resultFile)
		if err != nil {
			return err
		}
		defer outputfile.Close()

		r1, w1 := io.Pipe()
		r2, w2 := io.Pipe()

		sr1 := shell.NewShellRunner(nil, w1, nil, true)
		sr2 := shell.NewShellRunner(r1, w2, nil, true)
		sr3 := shell.NewShellRunner(r2, outputfile, nil, true)

		if err := sr1.RunShell("bcftools", "mpileup", "--seed", "42", "-Ou", "--max-depth", "100000", "-f", genomeFile, "-R", annotationFile, bamFile); err != nil {
			return err
		}
		if err := sr2.RunShell("bcftools", "call", "-Ou", "-mv"); err != nil {
			return err
		}
		if err := sr3.RunShell("bcftools", "filter", "-g3", "-G10", "-e", "QUAL<15 || (DP4[0]+DP4[1] == 0) || DP < 10"); err != nil {
			return err
		}

		if err := sr1.Wait(); err != nil {
			return err
		}
		w1.Close()

		if err := sr2.Wait(); err != nil {
			return err
		}
		w2.Close()

		if err := sr3.Wait(); err != nil {
			return err
		}
		return nil
	}
}

func runBedtoolsIntersectMpileupAndGff(annotationFile, pileupFile, resultFile string) func() error {
	return func() error {
		if err := pathForFile(resultFile); err != nil {
			return err
		}

		outputfile, err := os.Create(resultFile)
		if err != nil {
			return err
		}
		defer outputfile.Close()

		sr := shell.NewShellRunner(nil, outputfile, nil)
		return sr.RunShell("bedtools", "intersect", "-wa", "-wb", "-a", annotationFile, "-b", pileupFile)
	}
}

func generateResultTable(matureMirnaSeqs, hairpinMirnaSeqs, genomeFile, intersectTsvFile, resultFile string) func() error {
	return func() error {
		if err := pathForFile(resultFile); err != nil {
			return err
		}
		const (
			chrCol           = "CHROM"
			typeCol          = "TYPE"
			startCol         = "START"
			endCol           = "END"
			strandCol        = "STRAND"
			annCol           = "ANNOT"
			posCol           = "POS"
			refCol           = "REF"
			altCol           = "ALT"
			descCol          = "DESC"
			readsCountCol    = "READS"
			refReadsCountCol = "REF_READS"
			altReadsCountCol = "ALT_READS"
			vafCol           = "VAF"
			refSeqCol        = "REF_SEQ"
			refSeqPosCol     = "REF_SEQ_POS"
			altSeqCol        = "ALT_SEQ"
		)

		columns := []string{chrCol, typeCol, startCol, endCol, strandCol, annCol, posCol, refCol, altCol,
			readsCountCol, refReadsCountCol, altReadsCountCol, vafCol, refSeqPosCol, refSeqCol}

		columnsIndexMap := map[string]int{
			chrCol:    0,
			typeCol:   2,
			startCol:  3,
			endCol:    4,
			strandCol: 6,
			annCol:    8,
			posCol:    10,
			refCol:    12,
			altCol:    13,
			descCol:   16,
		}

		csvData, err := readCsv(intersectTsvFile, '\t')
		if err != nil {
			return err
		}

		rf, err := os.Create(resultFile)
		if err != nil {
			return err
		}
		defer rf.Close()

		mirnaSequences, err := seqsMapFromFile(matureMirnaSeqs)
		if err != nil {
			return err
		}

		hairpinSequences, err := seqsMapFromFile(hairpinMirnaSeqs)
		if err != nil {
			return err
		}

		if _, err := rf.WriteString(strings.Join(columns, "\t") + "\n"); err != nil {
			return err
		}

		linesMap := make(map[string]struct{})
		for _, line := range csvData {
			lineCol := func(colName string) string {
				return line[columnsIndexMap[colName]]
			}

			chrName := lineCol(chrCol)
			t := lineCol(typeCol)
			start := strToInt(lineCol(startCol))
			end := strToInt(lineCol(endCol))
			strand := lineCol(strandCol)
			mirnaName := strings.Split(strings.Split(lineCol(annCol), "Name=")[1], ";")[0]
			pos := strToInt(lineCol(posCol))
			ref := dnaToRna(lineCol(refCol))
			alt := dnaToRna(lineCol(altCol))

			readsCount, refReadsCount, altReadsCount, err := readsDescriptionFromLine(lineCol(descCol))
			if err != nil {
				return fmt.Errorf("reads don't match for miRNA %s in pos %d: %w", mirnaName, pos, err)
			}

			vaf := float64(altReadsCount) / float64(readsCount)

			var refSeq *linear.Seq
			switch t {
			case "miRNA":
				refSeq = mirnaSequences[mirnaName]
			case "miRNA_primary_transcript":
				refSeq = hairpinSequences[mirnaName]
			default:
				return fmt.Errorf("bad mirna sequence type: %s", t)
			}

			complRef := charToAlphabetLetter(ref)
			var seqPos int
			switch strand {
			case "+":
				seqPos = pos - start
			case "-":
				seqPos = end - pos
				complRef, err = complementRNA(complRef)
				if err != nil {
					return err
				}
			default:
				return fmt.Errorf("bad mirna strand type: %s", t)
			}

			if actualRef := refSeq.Seq[seqPos]; complRef != actualRef {
				return fmt.Errorf("bad ref seq for miRNA %s in pos %d: ref: %s, actual: %s", mirnaName, seqPos, string(complRef), string(actualRef))
			}

			resultLineContent := [][]byte{
				[]byte(chrName), []byte(t),
				intToByte(start), intToByte(end),
				[]byte(strand), []byte(mirnaName),
				intToByte(pos), []byte(ref),
				[]byte(alt), intToByte(readsCount),
				intToByte(refReadsCount), intToByte(altReadsCount),
				flaotToByte(vaf), intToByte(seqPos),
				[]byte(refSeq.String()),
			}
			resultLine := bytes.Join(resultLineContent, []byte("\t"))

			if _, f := linesMap[string(resultLine)]; f {
				continue
			}
			linesMap[string(resultLine)] = struct{}{}

			if _, err := rf.Write(resultLine); err != nil {
				return err
			}
			if _, err := rf.Write([]byte("\n")); err != nil {
				return err
			}
		}
		return nil
	}
}

func generateHeatmap(table, resultHeatMapImage string) func() error {
	return func() error {
		if err := pathForFile(resultHeatMapImage); err != nil {
			return err
		}
		sr := shell.NewShellRunner(nil, nil, nil)
		return sr.RunShell("python3", "scripts/heatmap.py", table, resultHeatMapImage)
	}
}

func generateBarplot(table, resultBarplotImage string) func() error {
	return func() error {
		if err := pathForFile(resultBarplotImage); err != nil {
			return err
		}
		sr := shell.NewShellRunner(nil, nil, nil)
		return sr.RunShell("python3", "scripts/barplot.py", table, resultBarplotImage)
	}
}

func readsDescriptionFromLine(l string) (int, int, int, error) {
	var readsCount, refReadsCount, altReadsCount int
	for _, colValue := range strings.Split(l, ";") {
		valueList := strings.Split(colValue, "=")
		col, v := valueList[0], valueList[1]
		switch col {
		case "DP":
			readsCount = strToInt(v)
		case "DP4":
			counts := strings.Split(v, ",")
			if len(counts) < 4 {
				return 0, 0, 0, fmt.Errorf("len of counts fields in DP4 is less then 4: '%s'", l)
			}

			refPrim, refSub := strToInt(counts[0]), strToInt(counts[1])
			altPrim, altSub := strToInt(counts[2]), strToInt(counts[3])

			if (refPrim != 0 && refSub != 0) || (altPrim != 0 && altSub != 0) || (refPrim != 0 && altSub != 0) || (altPrim != 0 && refSub != 0) {
				log.Printf("counts in DP4 are not on the same strands: '%s'\n", l)
			}

			refReadsCount = refPrim + refSub
			altReadsCount = altPrim + altSub
		}
	}
	if readsCount != (refReadsCount + altReadsCount) {
		return 0, 0, 0, fmt.Errorf("reads count don't match: full: %d, ref: %d, alt: %d", readsCount, refReadsCount, altReadsCount)
	}
	return readsCount, refReadsCount, altReadsCount, nil
}

func seqsMapFromFile(filePath string) (map[string]*linear.Seq, error) {
	f, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	sc := newRNAFastaReader(f)
	seqsMap := make(map[string]*linear.Seq)
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		seqsMap[s.ID] = s
	}
	if err := sc.Error(); err != nil {
		return nil, err
	}
	return seqsMap, nil
}

func chrMapper(tsvMapperFileName string) (map[string]string, error) {
	mapperData, err := readCsv(tsvMapperFileName, '\t')
	if err != nil {
		return nil, err
	}
	chrMapper := make(map[string]string, len(mapperData))
	for i, line := range mapperData {
		if len(line) < 2 {
			return nil, fmt.Errorf("malformed mapper on line %d", i+1)
		}
		chrMapper[line[0]] = line[1]
	}
	return chrMapper, nil
}
