package shell

import (
	"bufio"
	"io"
	"log"
	"os/exec"
	"strings"
)

type ShellRunner struct {
	stdin    io.Reader
	stdout   io.Writer
	stderr   io.Writer
	cmd      *exec.Cmd
	dontWait bool
}

func NewShellRunner(stdin io.Reader, stdout, sterr io.Writer, dontWait ...bool) *ShellRunner {
	var dw bool
	if len(dontWait) > 0 {
		dw = true
	}
	return &ShellRunner{
		stdin:    stdin,
		stdout:   stdout,
		stderr:   sterr,
		dontWait: dw,
	}
}

func (r *ShellRunner) RunShell(app string, args ...string) error {
	log.Printf("executing command '%s %s'\n", app, strings.Join(args, " "))
	r.cmd = exec.Command(app, args...)
	if r.stdin != nil {
		r.cmd.Stdin = r.stdin
	}
	var stdout, stderr io.ReadCloser
	if r.stdout != nil {
		r.cmd.Stdout = r.stdout
	} else {
		var err error
		stdout, err = r.cmd.StdoutPipe()
		if err != nil {
			return err
		}
	}

	if r.stderr != nil {
		r.cmd.Stderr = r.stderr
	} else {
		var err error
		stderr, err = r.cmd.StderrPipe()
		if err != nil {
			return err
		}

	}

	if err := r.cmd.Start(); err != nil {
		return err
	}

	go checkOutput(stdout)
	go checkOutput(stderr)

	if !r.dontWait {
		if err := r.Wait(); err != nil {
			return err
		}
	}
	return nil
}
func (r *ShellRunner) Wait() error {
	return r.cmd.Wait()
}

func checkOutput(reader io.ReadCloser) {
	if reader == nil {
		return
	}
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		log.Println(scanner.Text())
	}
}
