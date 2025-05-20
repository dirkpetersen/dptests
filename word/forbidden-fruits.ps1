# Script to highlight terms in Word documents as per CONVENTIONS.md

[CmdletBinding()]
param (
    [Parameter(Mandatory=$true, Position=0)]
    [string]$FolderPath,

    [Parameter(Mandatory=$true, Position=1)]
    [string]$TermsFilePath
)

# Load the Word interop assembly
Add-Type -AssemblyName Microsoft.Office.Interop.Word

# Function to load terms from a file
function Load-TermsFromFile {
    [CmdletBinding()]
    param (
        [Parameter(Mandatory=$true)]
        [string]$FilePath
    )
    try {
        if (Test-Path -Path $FilePath -PathType Leaf) {
            # Read all lines and filter out empty ones
            $terms = Get-Content -Path $FilePath -Raw | 
                     ForEach-Object { $_ -split "`r`n|\r|\n" } | 
                     Where-Object { $_.Trim() -ne "" }
            Write-Host "Loaded $(($terms | Measure-Object).Count) terms from file."
            return $terms
        } else {
            Write-Error "Terms file not found: $FilePath"
            return $null
        }
    }
    catch {
        Write-Error "Error reading terms file '$FilePath': $($_.Exception.Message)"
        return $null
    }
}

# Function to generate term variations
function Generate-TermVariations {
    [CmdletBinding()]
    param (
        [Parameter(Mandatory=$true)]
        [string]$Term
    )
    $variations = [System.Collections.ArrayList]::new()
    $trimmedTerm = $Term.Trim()

    if ($trimmedTerm.Length -gt 0) {
        [void]$variations.Add($trimmedTerm)
        if ($trimmedTerm.Contains(" ")) {
            # Create variations and add them separately
            $hyphenated = $trimmedTerm -replace " ", "-"
            $noSpaces = $trimmedTerm -replace " ", ""
            [void]$variations.Add($hyphenated)
            [void]$variations.Add($noSpaces)
        }
    }
    return $variations | Select-Object -Unique
}

# Function to process a single Word document
function Process-WordDocument {
    [CmdletBinding()]
    param (
        [Parameter(Mandatory=$true)]
        [string]$DocumentPath,

        [Parameter(Mandatory=$true)]
        [array]$TermsToHighlight,

        [Parameter(Mandatory=$true)]
        [System.Object]$WordApplication, # Word.Application COM object

        [Parameter(Mandatory=$true)]
        [int]$HighlightColor,
        
        [Parameter(Mandatory=$false)]
        [hashtable]$TermCountsRef
    )

    $doc = $null
    try {
        Write-Host "Opening document: $DocumentPath"
        $doc = $WordApplication.Documents.Open($DocumentPath)

        $contentRange = $doc.Content
        $find = $contentRange.Find

        foreach ($term in $TermsToHighlight) {
            if (-not ($term -and $term.Trim())) { continue } # Skip empty terms

            # Reset find parameters for each term
            $find.ClearFormatting()
            $find.Replacement.ClearFormatting() # Not strictly needed for find, but good practice

            $find.Text = $term
            $find.Forward = $true
            $find.Wrap = $wdFindStop # 0 = wdFindStop
            $find.Format = $false
            $find.MatchCase = $false
            $find.MatchWholeWord = $true
            $find.MatchWildcards = $false
            $find.MatchSoundsLike = $false
            $find.MatchAllWordForms = $false

            # Count occurrences for this term in this document
            $termCount = 0
            
            # Execute Find in a loop for all occurrences
            # The $contentRange is updated by $find.Execute() to the found range.
            # The next Execute() call continues from the end of the previously found range.
            while ($find.Execute()) {
                # Check if the found range ($contentRange) needs highlighting
                if ($contentRange.HighlightColorIndex -ne $HighlightColor) {
                    $contentRange.HighlightColorIndex = $HighlightColor
                    $termCount++
                }
                # $contentRange is now the found instance.
                # To prevent re-finding the same instance in some edge cases or if search doesn't advance,
                # one might collapse the range to its end: $contentRange.Collapse([Microsoft.Office.Interop.Word.WdCollapseDirection]::wdCollapseEnd)
                # However, standard Find.Execute() loop handles advancing correctly.
            }
            
            # Update the term counts in the reference hashtable
            if ($TermCountsRef -ne $null) {
                # Always update the hashtable, even if count is 0
                if ($TermCountsRef.ContainsKey($term)) {
                    $TermCountsRef[$term] += $termCount
                } else {
                    $TermCountsRef[$term] = $termCount
                }
            }
        }

        $doc.Save()
        Write-Host "Successfully processed and saved: $($doc.Name)"
    }
    catch {
        Write-Error "Error processing document '$($DocumentPath)': $($_.Exception.Message)"
        # Document might be partially highlighted. Save() is in try block.
        # If error occurs, changes up to that point might be saved if Save() was reached,
        # or not saved if error was before Save().
    }
    finally {
        if ($doc -ne $null) {
            $doc.Close()
            [System.Runtime.InteropServices.Marshal]::ReleaseComObject($doc) | Out-Null
            $doc = $null # Ensure it's nulled for GC
        }
    }
}

# --- Main Script Execution ---

# Define Word constants - using integer values directly
$wdYellow = 7 # wdYellow constant value
$wdFindStop = 0 # wdFindStop constant value

# 1. Load terms from the file
Write-Host "Loading terms from: $TermsFilePath"
$rawTerms = Load-TermsFromFile -FilePath $TermsFilePath
if (-not $rawTerms) {
    Write-Error "No terms loaded or terms file is empty. Exiting."
    exit 1
}

# 2. Generate all term variations and de-duplicate
$allTermsToHighlight = [System.Collections.ArrayList]::new()
foreach ($term in $rawTerms) {
    Write-Host "Processing term: $term"
    $variations = Generate-TermVariations -Term $term
    foreach ($variation in $variations) {
        [void]$allTermsToHighlight.Add($variation)
    }
}
$uniqueTermsToHighlight = $allTermsToHighlight | Select-Object -Unique | Where-Object { $_ -ne $null -and $_.Trim() -ne "" }

if ($uniqueTermsToHighlight.Count -eq 0) {
    Write-Warning "No valid terms to highlight after processing variations. Exiting."
    exit 0
}
Write-Host "Generated $(($uniqueTermsToHighlight | Measure-Object).Count) unique terms and variations to search for."
Write-Host "First 10 terms: $($uniqueTermsToHighlight | Select-Object -First 10 | ForEach-Object { "'$_'" } | Join-String -Separator ', ')"

# 3. Initialize Word Application
$wordApp = $null
try {
    Write-Host "Starting Microsoft Word application..."
    $wordApp = New-Object -ComObject Word.Application
    $wordApp.Visible = $false # Run Word in the background. Set to $true for debugging.
}
catch {
    Write-Error "Microsoft Word could not be started. Please ensure it is installed. Error: $($_.Exception.Message)"
    exit 1
}

# 4. Get Word documents from the specified folder
Write-Host "Searching for Word documents in: $FolderPath"
$wordFiles = Get-ChildItem -Path $FolderPath -Recurse -Include "*.doc", "*.docx" | Where-Object { $_.Name -notlike "~*" } # Exclude temporary Word files

if ($wordFiles.Count -eq 0) {
    Write-Warning "No Word documents found in '$FolderPath' or its subfolders."
} else {
    Write-Host "Found $($wordFiles.Count) Word document(s) to process."
    
    # Create a hashtable to track term counts across all documents
    $termCounts = @{}
    
    # 5. Process each document
    foreach ($fileInfo in $wordFiles) {
        Process-WordDocument -DocumentPath $fileInfo.FullName -TermsToHighlight $uniqueTermsToHighlight -WordApplication $wordApp -HighlightColor $wdYellow -TermCountsRef $termCounts
    }
    
    # 6. Display term counts in a table, sorted by frequency
    if ($termCounts.Count -gt 0) {
        Write-Host "`n--- Highlighted Terms Summary ---"
        
        # Make sure all terms from uniqueTermsToHighlight are in the counts hashtable
        foreach ($term in $uniqueTermsToHighlight) {
            if (-not $termCounts.ContainsKey($term)) {
                $termCounts[$term] = 0
            }
        }
        
        $termCounts.GetEnumerator() | 
            Sort-Object -Property Value -Descending | 
            Select-Object @{Name="Term"; Expression={$_.Key}}, @{Name="Count"; Expression={$_.Value}} | 
            Format-Table -AutoSize
    } else {
        Write-Host "No terms were highlighted in any documents."
    }
}

# 7. Clean up: Close Word application and release COM object
if ($wordApp -ne $null) {
    Write-Host "Closing Microsoft Word application."
    $wordApp.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wordApp) | Out-Null
    $wordApp = $null # Ensure it's nulled for GC
}

# Force garbage collection to help release COM objects, especially in interactive sessions
[System.GC]::Collect()
[System.GC]::WaitForPendingFinalizers()

Write-Host "Script finished."
