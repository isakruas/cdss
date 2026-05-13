# CDSS - Chirp-DSSS-Polar Digital Modem

CDSS is a Fortran 2018 implementation of a digital modem based on
chirp-assisted direct-sequence spread spectrum (DSSS), short-block Polar
coding, soft-decision reception, and iterative receiver refinement.

The project provides a command-line modem, reusable Fortran modules, waveform
file I/O, and a regression test suite intended to support maintenance, peer
review, and reproducible scientific experiments.

## Main Capabilities

- Encode binary payloads into complex baseband IQ waveforms.
- Decode synchronized IQ or audio captures back into payload bits.
- Synchronize captures by matched preamble correlation.
- Read and write `cf32`, mono WAV audio, and stereo WAV I/Q files.
- Run unit, integration, channel, and CLI regression tests with diagnostic logs.

## Modem Overview

The transmitter builds a protected frame, applies channel coding, spreads the
coded symbols, and synthesizes a chirp-assisted waveform. The receiver performs
preamble acquisition, soft demodulation, reliability weighting, Polar decoding,
and CRC validation.

```text
Payload bits
  -> CRC-16 protected header
  -> Polar channel coding
  -> BPSK symbol mapping
  -> Gold-code DSSS spreading
  -> chirp carrier synthesis
  -> complex baseband waveform
```

```text
Received waveform
  -> preamble synchronization
  -> chip-level soft demodulation
  -> uncertainty/reliability weighting
  -> Polar SCL decoding
  -> CRC validation
  -> recovered payload bits
```

The frame layout is:

```text
31-bit preamble | Polar-coded 48-bit header | Polar-coded payload
```

The header stores the payload length and a CRC-16/CCITT-FALSE checksum of the
payload bits.

## Repository Layout

```text
cdss/
  app/
    cdss_cli.f90              command-line interface entry point
  src/
    api/                      public modem API
    channel/                  AWGN, CFO, drift, tones, and clock error models
    codec/                    Polar encoder and SCL decoder
    core/                     constants, kinds, and shared derived types
    receiver/                 synchronization, demodulation, and refinement
    util/                     bit utilities, CRC, RNG, FFTW, and file I/O
    waveform/                 PN sequences, chirp carrier, and synthesis
  test/                       automated test programs
  fpm.toml                    Fortran Package Manager project definition
  Makefile                    documented build, test, install, and run targets
  LICENSE                     project license
```

## Requirements

- `gfortran` with Fortran 2018 support
- `fpm` (Fortran Package Manager)
- FFTW3 development files (`fftw3.f03` and `libfftw3`)
- `make`
- POSIX shell utilities used by the Makefile and CLI tests

On Debian or Ubuntu systems:

```bash
sudo apt install gfortran fpm libfftw3-dev make
```

## Build

The Makefile exports the compiler and linker flags required by the FFTW Fortran
interface.

```bash
make build       # optimized release build
make debug       # debug profile build
make install     # install cdss to PREFIX/bin, default ~/.local/bin
make clean       # remove FPM build artifacts
```

The installation prefix can be changed with:

```bash
make install PREFIX=/opt/cdss
```

## Command-Line Usage

Print command help:

```bash
cdss help
```

Print modem geometry and default parameters as JSON:

```bash
cdss info
```

### Modulate

Encode a bit string into a waveform file:

```bash
cdss mod --bit-string 101100111000 --output frame.cf32 --output-format cf32
```

Supported output formats:

- `cf32`: interleaved 32-bit floating-point I/Q samples
- `wav-audio`: mono 16-bit PCM audio projection
- `wav-iq`: stereo 16-bit PCM I/Q, with I on the left channel and Q on the right

### Synchronize

Detect the preamble in a capture and write a frame-aligned waveform:

```bash
cdss sync --input capture.cf32 --output synced.cf32 --bits 64
```

### Demodulate

Decode a frame-aligned waveform:

```bash
cdss demod --input synced.cf32 --input-format cf32 --max-payload-bits 128
```

## Default Parameters

Default modem constants are defined in `src/core/cdss_types.f90` and
`src/core/cdss_constants.f90`.

| Parameter | Default |
|---|---:|
| Sample rate | 6000 Hz |
| Carrier frequency | 1500 Hz |
| Chip rate | 250 chips/s |
| Spreading factor | 64 chips/coded bit |
| Coded-bit rate | 3.90625 bit/s |
| Chirp slope | 50 Hz/s |
| Polar code | N=128, K=64 |
| SCL list size | 8 |
| Preamble length | 31 bits |
| Header length | 48 bits |
| Refinement iterations | 2 |

## Tests

Run the complete regression suite with:

```bash
make test
```

The tests print diagnostic values for maintainers and reviewers, including CRC
values, decoded lengths, synchronization scores, LLR ranges, channel model
measurements, waveform round-trip error, CLI JSON output, and expected failure
exit statuses.

The suite covers:

- bit conversion, bit-string parsing, CRC-16, and frame header handling
- Polar encoding and SCL decoding
- PN and Gold-code sequence generation
- waveform synthesis, carrier generation, envelopes, and deterministic output
- AWGN, CFO, drift, clock error, guard padding, and tone interference
- `cf32`, `wav-audio`, and `wav-iq` I/O round trips
- synchronization, soft demodulation, uncertainty weighting, and refinement
- zero-length, boundary-length, and multi-block payloads
- high-SNR integration and deterministic channel regressions
- CLI success paths and expected CLI error paths

## Development Workflow

Typical local workflow:

```bash
make debug
make test
make run ARGS='info'
```

## License

CDSS is distributed under the GNU General Public License v3.0 or later. See
`LICENSE` for the full license text.
