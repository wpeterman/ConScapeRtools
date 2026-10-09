param(
    [Parameter(Mandatory=$true)][string]$Julia,
    [Parameter(Mandatory=$true)][string]$Repository,
    [Parameter(Mandatory=$true)][string]$OutputDirectory
)
# Run alone: other benchmark/check jobs compete for RAM and CPU.
# PeakWorkingSet64 is the OS high-water mark sampled every 100 ms. It includes
# loading, input construction and the small warm-up, for one Julia process.
$ErrorActionPreference = 'Stop'
$benchRoot = (Resolve-Path -LiteralPath $Repository).Path
if (Test-Path -LiteralPath $OutputDirectory) { throw 'Choose a new output directory.' }
New-Item -ItemType Directory -Path $OutputDirectory | Out-Null
$cases = @()
foreach ($shape in @(@(50, 1), @(80, 1), @(120, 0.02))) {
    foreach ($method in @('separate', 'shared', 'streamed')) {
        $cases += [pscustomobject]@{Side=$shape[0]; Density=$shape[1]; Method=$method; Blas=1; Block=128}
    }
}
foreach ($block in @(16, 64, 256)) {
    $cases += [pscustomobject]@{Side=80; Density=1; Method='streamed'; Blas=1; Block=$block}
}
$cases += [pscustomobject]@{Side=80; Density=1; Method='streamed'; Blas=2; Block=128}
$peaks = @()
foreach ($case in $cases) {
    $label = "$($case.Method)-$($case.Side)-$($case.Density)-blas$($case.Blas)-block$($case.Block)"
    Write-Output "Running $label"
    $out = Join-Path $OutputDirectory "$label.csv"
    $err = Join-Path $OutputDirectory "$label.stderr"
    $driver = Join-Path $benchRoot 'inst/examples/benchmark_julia.jl'
    $argsList = @('--startup-file=no', ('"' + $driver + '"'), ('"' + $benchRoot + '"'),
        $case.Side, $case.Density, $case.Method, $case.Blas, 3, $case.Block)
    $proc = Start-Process -FilePath $Julia -ArgumentList $argsList -PassThru -WindowStyle Hidden `
        -RedirectStandardOutput $out -RedirectStandardError $err
    $peak = 0L
    while (-not $proc.HasExited) {
        $proc.Refresh()
        $peak = [Math]::Max($peak, $proc.PeakWorkingSet64)
        Start-Sleep -Milliseconds 100
    }
    $proc.WaitForExit()
    if ($proc.ExitCode -ne 0) { throw "Failed $label. See $err" }
    $peaks += [pscustomobject]@{Case=$label; PeakWorkingSetMB=$peak / 1e6}
    $peaks | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'peaks.csv') -NoTypeInformation
    Get-Content -LiteralPath $out | Select-Object -Last 3
    Write-Output "Peak working set: $([Math]::Round($peak / 1e6, 2)) MB"
}
