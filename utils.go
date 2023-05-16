package main

import (
	"archive/zip"
	"bytes"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"compress/gzip"

	"github.com/alex123012/final-thesis/pkg/shell"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/schollz/progressbar/v3"
)

func pathForFile(fileName string) error {
	filePath := filepath.Dir(fileName)
	if filePath == "" {
		return nil
	}
	if err := mkdir(filePath); err != nil {
		return err
	}
	return nil
}

func mkdir(dir string) error {
	if err := os.MkdirAll(dir, 0755); err != nil {
		return err
	}
	return nil
}

func fileExists(fileNames ...string) bool {
	for _, fileName := range fileNames {
		if _, err := os.Stat(fileName); errors.Is(err, os.ErrNotExist) {
			return false
		}
	}
	return true
}

func downloadMirnaAnnotations(filePath string) error {
	fileUrl := "https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3"
	return downloadFile(fileUrl, filePath)
}

func downloadMirnaSeqs(t string, filepath string) error {
	fileUrl := fmt.Sprintf("https://www.mirbase.org/ftp/CURRENT/%s.fa.gz", t)
	archive := fmt.Sprintf("%s.fa.gz", t)
	fileInArchive := fmt.Sprintf("%s.fa", t)

	if err := downloadFile(fileUrl, archive); err != nil {
		return err
	}
	if err := extractFileFromGz(archive, fileInArchive, filepath); err != nil {
		return err
	}
	if err := os.Remove(archive); err != nil {
		return err
	}
	return nil
}

func downloadGenomeSequence(genomeVersion, genomeFile string) (string, error) {
	zipFile := genomeFile + ".zip"
	if fileExists(zipFile) {
		return zipFile, nil
	}

	fullURLFile := fmt.Sprintf("https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/%s/download?include_annotation_type=GENOME_GFF&filename=%s", genomeVersion, zipFile)

	return zipFile, downloadFile(fullURLFile, zipFile)
}

func downloadFile(fileUrl, filePathToSave string) error {
	req, err := http.NewRequest("GET", fileUrl, nil)
	if err != nil {
		return err
	}
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return err
	}

	defer resp.Body.Close()

	bar := progressbar.DefaultBytes(
		resp.ContentLength,
		fmt.Sprintf("downloading %s", filePathToSave),
	)

	file, err := os.Create(filePathToSave)
	if err != nil {
		return err
	}
	defer file.Close()

	_, err = io.Copy(io.MultiWriter(file, bar), resp.Body)
	return err
}

func extractFileFromZip(archiveName, fileInArchive, resultFileName string) error {
	r, err := zip.OpenReader(archiveName)
	if err != nil {
		return err
	}
	defer r.Close()

	rc, err := r.Open(fileInArchive)
	if err != nil {
		return err
	}
	defer rc.Close()

	rf, err := os.Create(resultFileName)
	if err != nil {
		return err
	}
	defer rf.Close()

	if _, err := io.Copy(rf, rc); err != nil {
		return err
	}
	return nil
}

func extractFileFromGz(archiveName, fileInArchive, resultFileName string) error {
	f, err := os.Open(archiveName)
	if err != nil {
		return err
	}
	defer f.Close()

	rf, err := os.Create(resultFileName)
	if err != nil {
		return err
	}
	defer rf.Close()

	gr, err := gzip.NewReader(f)
	if err != nil {
		return err
	}
	if _, err := io.Copy(rf, gr); err != nil {
		return err
	}
	return nil
}

func readCsv(fileName string, sep rune) ([][]string, error) {
	csvFile, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer csvFile.Close()

	reader := csv.NewReader(csvFile)
	reader.Comma = sep
	reader.FieldsPerRecord = -1

	csvData, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}
	return csvData, nil
}

func newDNAFastaReader(f io.Reader) *seqio.Scanner {
	return newFastaReader(f, alphabet.DNAredundant)
}

func newRNAFastaReader(f io.Reader) *seqio.Scanner {
	return newFastaReader(f, alphabet.RNAredundant)
}

func newFastaReader(f io.Reader, a alphabet.Alphabet) *seqio.Scanner {
	return seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, a)))
}

func intToByte(v int) []byte {
	return []byte(strconv.Itoa(v))
}

func flaotToByte(v float64) []byte {
	return []byte(fmt.Sprintf("%f", v))
}

func strToInt(v string) int {
	i, _ := strconv.Atoi(v)
	return i
}

func complementRNA(n alphabet.Letter) (alphabet.Letter, error) {
	c, ok := alphabet.RNAredundant.Complement(alphabet.Letter(n))
	if !ok {
		return 0, fmt.Errorf("not a valid DNA nucleotide: %b", n)
	}
	return c, nil
}

func dnaToRna(ns string) string {
	splitted := strings.Split(ns, ",")
	res := make([]string, 0, len(splitted))
	for _, n := range splitted {
		if n == "T" {
			n = "U"
		}
		res = append(res, n)
	}
	return strings.Join(res, ",")
}

func charToAlphabetLetter(c string) alphabet.Letter {
	return alphabet.Letter([]byte(c)[0])
}

func seqToBytes(s *linear.Seq) []byte {
	b := bytes.NewBuffer(nil)
	b.WriteString(fasta.DefaultIDPrefix)
	b.WriteString(s.ID)
	b.WriteByte(' ')
	b.WriteString(s.Desc)
	b.WriteByte('\n')
	b.Write(alphabet.LettersToBytes(s.Seq))
	b.WriteByte('\n')
	return b.Bytes()
}

func gzipFile(fileName string) func() error {
	return func() error {
		sr := shell.NewShellRunner(nil, nil, nil)
		return sr.RunShell("gzip", fileName)
	}
}
