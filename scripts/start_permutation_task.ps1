# Start (or resume) the Stage 4 permutation run as a one-off Windows scheduled task,
# independent of any terminal or Claude session. The task runs
# scripts/run_permutation_task.cmd, which deletes the task when the runner exits.
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$name = "TCGA_BRCA_permutation"

if (Get-ScheduledTask -TaskName $name -ErrorAction SilentlyContinue) {
    Write-Output "Task $name already exists (state: $((Get-ScheduledTask -TaskName $name).State)). Not starting a second one."
    exit 0
}
$cmd = Join-Path $root "scripts\run_permutation_task.cmd"
$action = New-ScheduledTaskAction -Execute "conhost.exe" -Argument "--headless cmd.exe /c `"`"$cmd`"`"" -WorkingDirectory $root
$trigger = New-ScheduledTaskTrigger -Once -At (Get-Date).AddMinutes(1)
$settings = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit (New-TimeSpan -Hours 24) -MultipleInstances IgnoreNew
$principal = New-ScheduledTaskPrincipal -UserId "$env:USERDOMAIN\$env:USERNAME" -LogonType Interactive -RunLevel Limited
Register-ScheduledTask -TaskName $name -Action $action -Trigger $trigger -Settings $settings `
    -Principal $principal -Description "One-off Stage 4 permutation run; deletes itself when done." | Out-Null
Start-ScheduledTask -TaskName $name
Write-Output "Started task $name"
