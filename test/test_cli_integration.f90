!> @file test_cli_integration.f90
!! @brief Integration tests for the installed command-line executable.
program test_cli_integration
  implicit none

  character(len=*), parameter :: info_path = "/tmp/cdss_cli_info.json"
  character(len=*), parameter :: mod_path = "/tmp/cdss_cli_frame.cf32"
  character(len=*), parameter :: sync_path = "/tmp/cdss_cli_sync.cf32"
  character(len=*), parameter :: demod_path = "/tmp/cdss_cli_demod.json"
  character(len=*), parameter :: simulate_path = "/tmp/cdss_cli_simulate.json"
  character(len=*), parameter :: usage_path = "/tmp/cdss_cli_usage.txt"
  character(len=*), parameter :: bad_path = "/tmp/cdss_cli_bad.txt"
  character(len=4096) :: cmd
  character(len=:), allocatable :: text
  integer :: stat

  print '(a)', "test_cli_integration: starting CLI checks"

  call run_cli("> " // usage_path // " 2>&1", stat)
  print '(a,i0)', "  no-argument usage exit status=", stat
  if (stat /= 2) error stop 1
  call read_text(usage_path, text)
  call require_contains(text, 'usage: cdss <command> [options]', "usage header")
  call require_contains(text, 'sync --input <path>', "usage sync command")
  call require_contains(text, '--sync-search', "usage sync-search option")
  call require_contains(text, '--tone-frequency-hz X', "usage tone option")
  call require_not_contains(text, 'STOP', "usage should not print STOP")

  call run_cli("info > " // info_path, stat)
  if (stat /= 0) error stop 2
  call read_text(info_path, text)
  call require_contains(text, '"sample_rate_hz"', "info sample_rate_hz")
  call require_contains(text, '"spreading_factor"', "info spreading_factor")
  print '(a,a)', "  info output=", trim(text)

  call run_cli("mod --bit-string 101100111000 --output " // mod_path // &
       " --output-format cf32 > /tmp/cdss_cli_mod.json", stat)
  if (stat /= 0) error stop 3

  call run_cli("sync --input " // mod_path // " --output " // sync_path // &
       " --bits 12 > /tmp/cdss_cli_sync.json", stat)
  if (stat /= 0) error stop 4
  call read_text("/tmp/cdss_cli_sync.json", text)
  call require_contains(text, '"detected_sample"', "sync detected_sample")
  call require_contains(text, '"frame_samples"', "sync frame_samples")
  print '(a,a)', "  sync output=", trim(text)

  call run_cli("demod --input " // mod_path // " --max-payload-bits 64 > " // demod_path, stat)
  if (stat /= 0) error stop 5
  call read_text(demod_path, text)
  call require_contains(text, '"crc_ok":true', "demod crc_ok")
  call require_contains(text, '"hard_bits":"101100111000"', "demod hard_bits")
  print '(a,a)', "  demod output=", trim(text)

  call run_cli("simulate --bits 16 --snr-db 45 --seed 7 --cfo-hz 0.001 " // &
       "--drift-hz-s 0.0001 --clock-ppm 0 --prefix-pad-s 0.01 " // &
       "--suffix-pad-s 0.005 --sync-search --tone-frequency-hz 1525 --tone-sir-db 35 > " // &
       simulate_path, stat)
  if (stat /= 0) error stop 6
  call read_text(simulate_path, text)
  call require_contains(text, '"crc_ok":true', "simulate crc_ok")
  call require_contains(text, '"bit_errors":0', "simulate bit_errors")
  call require_contains(text, '"tone_count":1', "simulate tone_count")
  call require_contains(text, '"sync_search":true', "simulate sync_search")
  call require_contains(text, '"sync_ok":true', "simulate sync_ok")
  call require_contains(text, '"prefix_pad_s":1.00000000E-02', "simulate prefix_pad_s")
  print '(a,a)', "  simulate output=", trim(text)

  call run_cli("unknown > " // bad_path // " 2>&1", stat)
  print '(a,i0)', "  unknown command exit status=", stat
  if (stat == 0) then
    print '(a)', "test_cli_integration: unknown command should fail"
    error stop 7
  end if

  call run_cli("mod --bit-string 101 --output /tmp/cdss_bad.out " // &
       "--output-format unsupported > " // bad_path // " 2>&1", stat)
  print '(a,i0)', "  unsupported format exit status=", stat
  if (stat == 0) then
    print '(a)', "test_cli_integration: unsupported output format should fail"
    error stop 8
  end if

  call run_cli("demod --input /tmp/cdss_missing_input.cf32 > " // bad_path // " 2>&1", stat)
  print '(a,i0)', "  missing input exit status=", stat
  if (stat == 0) then
    print '(a)', "test_cli_integration: missing input file should fail"
    error stop 9
  end if

  call delete_file(usage_path)
  call delete_file(info_path)
  call delete_file(mod_path)
  call delete_file(sync_path)
  call delete_file(demod_path)
  call delete_file(simulate_path)
  call delete_file(bad_path)
  call delete_file("/tmp/cdss_cli_mod.json")
  call delete_file("/tmp/cdss_cli_sync.json")
  print '(a)', "test_cli_integration: passed"

contains

  subroutine run_cli(args, exit_status)
    character(len=*), intent(in) :: args
    integer, intent(out) :: exit_status
    character(len=4096) :: command
    command = "sh -c 'bin=$(find build -path ""*/app/cdss"" -executable -type f " // &
         "-printf ""%T@ %p\n"" | sort -nr | head -1 | cut -d"" "" -f2-); " // &
         "test -n ""$bin"" && ""$bin"" " // args // "'"
    call execute_command_line(trim(command), exitstat=exit_status)
  end subroutine run_cli

  subroutine require_contains(haystack, needle, label)
    character(len=*), intent(in) :: haystack, needle, label
    if (index(haystack, needle) <= 0) then
      print '(a,a)', "test_cli_integration: missing expected token: ", label
      print '(a,a)', "  expected token=", needle
      print '(a,a)', "  output=", trim(haystack)
      error stop 10
    end if
  end subroutine require_contains

  subroutine require_not_contains(haystack, needle, label)
    character(len=*), intent(in) :: haystack, needle, label
    if (index(haystack, needle) > 0) then
      print '(a,a)', "test_cli_integration: unexpected token: ", label
      print '(a,a)', "  unexpected token=", needle
      print '(a,a)', "  output=", trim(haystack)
      error stop 12
    end if
  end subroutine require_not_contains

  subroutine read_text(path, content)
    character(len=*), intent(in) :: path
    character(len=:), allocatable, intent(out) :: content
    character(len=8192) :: buffer
    integer :: unit, ios, used
    buffer = ""
    used = 0
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(a,a)', "test_cli_integration: cannot open ", path
      error stop 11
    end if
    do
      read(unit, '(a)', iostat=ios) buffer(used + 1:)
      if (ios /= 0) exit
      used = len_trim(buffer)
      if (used >= len(buffer) - 1) exit
    end do
    close(unit)
    allocate(character(len=used) :: content)
    content = buffer(:used)
  end subroutine read_text

  subroutine delete_file(path)
    character(len=*), intent(in) :: path
    integer :: unit, ios
    open(newunit=unit, file=path, status='old', iostat=ios)
    if (ios == 0) close(unit, status='delete')
  end subroutine delete_file
end program test_cli_integration
