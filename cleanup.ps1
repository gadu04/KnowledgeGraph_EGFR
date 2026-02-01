# Cleanup old files after refactoring

Write-Host "`n=== 🧹 CLEANING UP OLD FILES ===`n" -ForegroundColor Cyan

# Move old Python files to archive
$oldFiles = @("BuildKG.py", "data.py", "Eval.py", "pipe.py")

foreach ($file in $oldFiles) {
    if (Test-Path $file) {
        $newName = [System.IO.Path]::GetFileNameWithoutExtension($file) + "_old.py"
        Move-Item $file "archive\$newName" -Force
        Write-Host "✅ Moved $file → archive\$newName" -ForegroundColor Green
    }
}

# Remove empty directories
$emptyDirs = @("backup", "output", "Data")

foreach ($dir in $emptyDirs) {
    if (Test-Path $dir) {
        $items = Get-ChildItem $dir -Recurse
        if ($items.Count -eq 0) {
            Remove-Item $dir -Force -Recurse
            Write-Host "✅ Removed empty directory: $dir" -ForegroundColor Green
        } else {
            Write-Host "⚠️  Directory $dir is not empty, skipping..." -ForegroundColor Yellow
        }
    }
}

Write-Host "`n=== ✅ CLEANUP COMPLETED ===`n" -ForegroundColor Green

# List remaining files in root
Write-Host "=== FILES IN ROOT DIRECTORY ===" -ForegroundColor Cyan
Get-ChildItem -File | Where-Object { $_.Extension -in @('.py', '.ipynb', '.yml', '.md') } | Select-Object Name | Format-Table -HideTableHeaders

Write-Host "`n=== NEXT STEPS ===" -ForegroundColor Yellow
Write-Host "1. Test the new structure: python scripts/build_kg.py"
Write-Host "2. Update notebook imports if needed"
Write-Host "3. Commit changes: git add . && git commit -m 'refactor: Clean code structure'"
Write-Host "4. Push to GitHub: git push origin main`n"
