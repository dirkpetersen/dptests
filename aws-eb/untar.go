package main

import (
    "archive/tar"
    "compress/gzip"
    "flag"
    "fmt"
    "io"
    "os"
    "path/filepath"
    "runtime"
    "strings"
    "sync"
)

func main() {
    // Parse command line arguments
    flag.Parse()
    args := flag.Args()
    if len(args) < 1 {
        fmt.Println("Usage: untar <directory> [parallel]")
        os.Exit(1)
    }
    rootDir := args[0]
    parallel := runtime.NumCPU() * 10
    if len(args) > 1 {
        fmt.Sscanf(args[1], "%d", &parallel)
    }

    // Channel to process files
    files := make(chan string, parallel)

    // WaitGroup for goroutines
    var wg sync.WaitGroup

    // Start workers
    for i := 0; i < parallel; i++ {
        wg.Add(1)
        go func() {
            defer wg.Done()
            for f := range files {
                untarFile(f)
            }
        }()
    }

    // Walk through the directory
    err := filepath.Walk(rootDir, func(path string, info os.FileInfo, err error) error {
        if err != nil {
            fmt.Printf("Error accessing path %q: %v\n", path, err)
            return nil // Continue walking despite the error
        }
        if info.IsDir() {
            return nil // Continue walking
        }
        if strings.HasSuffix(path, ".tar.gz") {
            files <- path
            return filepath.SkipDir // Skip deeper exploration of this directory
        }
        return nil
    })

    if err != nil {
        fmt.Printf("Error walking the path %q: %v\n", rootDir, err)
    }

    // Close channel and wait for goroutines to finish
    close(files)
    wg.Wait()
}

func untarFile(filePath string) {
    fmt.Println("Untarring:", filePath)
    file, err := os.Open(filePath)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer file.Close()

    gzipReader, err := gzip.NewReader(file)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer gzipReader.Close()

    tarReader := tar.NewReader(gzipReader)

    // Get the directory of the file
    dir := filepath.Dir(filePath)

    for {
        header, err := tarReader.Next()
        if err == io.EOF {
            break
        }
        if err != nil {
            fmt.Println(err)
            return
        }

        // Construct the full path for the file to be created
        fullPath := filepath.Join(dir, header.Name)

        switch header.Typeflag {
        case tar.TypeDir:
            // Create directory
            if err := os.MkdirAll(fullPath, 0755); err != nil {
                fmt.Println(err)
            }
        case tar.TypeReg:
            // Create file
            outFile, err := os.Create(fullPath)
            if err != nil {
                fmt.Println(err)
                continue
            }
            if _, err := io.Copy(outFile, tarReader); err != nil {
                fmt.Println(err)
            }
            outFile.Close()
        }
    }
}