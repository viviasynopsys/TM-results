$ok = 0; $bad = 0
Get-Content RESTORE_EXTENSIONS.md -Encoding UTF8 | ForEach-Object {
  if ($_ -match '^\| `([^`]+)` \| `([^`]+)` \| `[^`]*` \| `([0-9a-f]{64})` \|') {
    $src = $Matches[1]; $dst = $Matches[2]; $want = $Matches[3]
    if (Test-Path -LiteralPath $src) { Move-Item -LiteralPath $src -Destination $dst -Force }
    if (-not (Test-Path -LiteralPath $dst)) { Write-Output "MISSING $dst"; $bad++; return }
    $got = (Get-FileHash -LiteralPath $dst -Algorithm SHA256).Hash.ToLower()
    if ($got -eq $want) { $ok++; Write-Output "OK  $dst" } else { $bad++; Write-Output "BAD $dst" }
  }
}
Write-Output "restored $ok OK, $bad BAD"
