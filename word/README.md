# Forbidden Fruits 

Forbidden Fruit  is a Powershell script that reviews a folder of Word documents and drives MS Word to highlight all terms yellow that are kept in an external text file (one term per line) 


```
powershell ./forbidden-fruits.ps1 -FolderPath ./docs -TermsFilePath ./forbidden-fruits-terms.txt


```

there is also a Python implementation, however the docx package it uses is less reliable 

