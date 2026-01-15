#!/usr/bin/env python3
"""
fqview - A colorful FASTQ file viewer for the terminal

Displays FASTQ records with:
- Color-coded header fields (colon-separated parts)
- Sequence with optional coloring
- Quality scores as colored blocks
"""

import argparse
import sys
import gzip
from pathlib import Path


# ANSI color codes
class Colors:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    
    # Header field colors (rotating through these for colon-separated parts)
    HEADER_COLORS = [
        '\033[38;5;39m',   # Bright blue
        '\033[38;5;208m',  # Orange
        '\033[38;5;170m',  # Pink/magenta
        '\033[38;5;114m',  # Light green
        '\033[38;5;220m',  # Yellow
        '\033[38;5;147m',  # Light purple
        '\033[38;5;87m',   # Cyan
        '\033[38;5;203m',  # Coral
    ]
    
    # Sequence base colors
    BASE_COLORS = {
        'A': '\033[38;5;46m',   # Green
        'T': '\033[38;5;196m',  # Red
        'G': '\033[38;5;226m',  # Yellow
        'C': '\033[38;5;33m',   # Blue
        'N': '\033[38;5;240m',  # Gray
    }


def get_quality_color(phred_score: int) -> str:
    """
    Get a foreground color for a Phred quality score.
    Maps Phred 0-41+ to a color gradient from red (bad) to green (good).
    """
    # Quality thresholds and corresponding 256-color codes
    # Red -> Orange -> Yellow -> Light Green -> Green
    if phred_score < 10:
        # Very low quality: shades of red
        return f'\033[38;5;{196 + min(phred_score, 5)}m'  # Red tones
    elif phred_score < 20:
        # Low quality: orange to yellow
        colors = [208, 214, 220, 226, 227, 228, 229, 230, 190, 191]
        idx = min(phred_score - 10, len(colors) - 1)
        return f'\033[38;5;{colors[idx]}m'
    elif phred_score < 30:
        # Medium quality: yellow-green
        colors = [192, 149, 150, 151, 152, 114, 115, 116, 84, 85]
        idx = min(phred_score - 20, len(colors) - 1)
        return f'\033[38;5;{colors[idx]}m'
    else:
        # High quality: greens
        if phred_score < 35:
            return '\033[38;5;82m'   # Bright green
        elif phred_score < 40:
            return '\033[38;5;46m'   # Pure green
        else:
            return '\033[38;5;48m'   # Cyan-green (excellent)


def format_header(line: str, prefix_char: str = '@') -> str:
    """Format a FASTQ header line with colors for each colon-separated field."""
    if not line:
        return line
    
    # Remove the prefix character
    content = line[1:] if line.startswith(prefix_char) else line
    
    # Split by spaces first (instrument info vs description)
    parts = content.split(' ', 1)
    main_part = parts[0]
    description = parts[1] if len(parts) > 1 else None
    
    # Color the colon-separated fields in the main part
    fields = main_part.split(':')
    colored_fields = []
    for i, field in enumerate(fields):
        color = Colors.HEADER_COLORS[i % len(Colors.HEADER_COLORS)]
        colored_fields.append(f"{color}{field}{Colors.RESET}")
    
    result = f"{Colors.DIM}{prefix_char}{Colors.RESET}" + ':'.join(colored_fields)
    
    # Add description if present
    if description:
        # Also color colon-separated parts in description
        desc_fields = description.split(':')
        colored_desc = []
        for i, field in enumerate(desc_fields):
            color = Colors.HEADER_COLORS[(i + len(fields)) % len(Colors.HEADER_COLORS)]
            colored_desc.append(f"{color}{field}{Colors.RESET}")
        result += f" {':'.join(colored_desc)}"
    
    return result


def format_sequence(seq: str, color_bases: bool = True) -> str:
    """Format a sequence line, optionally with colored bases."""
    if not color_bases:
        return seq
    
    result = []
    for base in seq:
        color = Colors.BASE_COLORS.get(base.upper(), Colors.RESET)
        result.append(f"{color}{base}{Colors.RESET}")
    return ''.join(result)


def format_quality(qual: str) -> str:
    """Format quality scores as colored blocks."""
    result = []
    for char in qual:
        # Convert ASCII to Phred score (Phred+33 encoding, most common)
        phred = ord(char) - 33
        color = get_quality_color(phred)
        # Use a solid block character with foreground color
        result.append(f"{color}█{Colors.RESET}")
    return ''.join(result)


def print_quality_legend():
    """Print a legend explaining the quality color scale."""
    print(f"\n{Colors.BOLD}Quality Score Legend:{Colors.RESET}")
    print("Phred: ", end='')
    
    # Show sample scores with their colors
    for score in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
        color = get_quality_color(score)
        print(f"{color}██{Colors.RESET}{score:<2}", end='')
    
    print(f"\n       {'Low':<12}{'Medium':<12}{'High':<12}")
    print()


def open_fastq(filepath: str):
    """Open a FASTQ file, handling gzip compression if needed."""
    path = Path(filepath)
    if path.suffix == '.gz' or filepath.endswith('.fastq.gz') or filepath.endswith('.fq.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def read_fastq_records(filepath: str, num_records: int = None):
    """Generator that yields FASTQ records as tuples of (header, seq, plus, qual)."""
    count = 0
    with open_fastq(filepath) as f:
        while True:
            if num_records is not None and count >= num_records:
                break
            
            header = f.readline().rstrip('\n')
            if not header:
                break
                
            seq = f.readline().rstrip('\n')
            plus = f.readline().rstrip('\n')
            qual = f.readline().rstrip('\n')
            
            if not all([seq, plus, qual]):
                break
                
            yield (header, seq, plus, qual)
            count += 1


def main():
    parser = argparse.ArgumentParser(
        description='Display FASTQ files with colorful quality visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fqview reads.fastq -n 5          # Show first 5 records
  fqview reads.fastq.gz -n 10      # Works with gzipped files
  fqview reads.fq --no-seq-color   # Disable sequence coloring
  fqview reads.fastq --legend      # Show quality color legend
        """
    )
    
    parser.add_argument('fastq', help='Input FASTQ file (supports .gz)')
    parser.add_argument('-n', '--num-records', type=int, default=4,
                        help='Number of records to display (default: 4)')
    parser.add_argument('--no-seq-color', action='store_true',
                        help='Disable coloring of sequence bases')
    parser.add_argument('--legend', action='store_true',
                        help='Show quality score color legend')
    parser.add_argument('--raw-quality', action='store_true',
                        help='Show raw quality string alongside colored blocks')
    parser.add_argument('--wrap', type=int, default=0,
                        help='Wrap long sequences at this width (0=no wrap)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.fastq).exists():
        print(f"Error: File not found: {args.fastq}", file=sys.stderr)
        sys.exit(1)
    
    # Print legend if requested
    if args.legend:
        print_quality_legend()
    
    # Print header
    print(f"{Colors.BOLD}{'─' * 60}{Colors.RESET}")
    print(f"{Colors.BOLD}FASTQ Viewer: {args.fastq}{Colors.RESET}")
    print(f"{Colors.BOLD}{'─' * 60}{Colors.RESET}\n")
    
    # Process records
    try:
        for i, (header, seq, plus, qual) in enumerate(read_fastq_records(args.fastq, args.num_records)):
            print(f"{Colors.DIM}Record {i+1}:{Colors.RESET}")
            
            # Header line
            print(f"  {format_header(header, '@')}")
            
            # Sequence
            if args.wrap > 0:
                for j in range(0, len(seq), args.wrap):
                    chunk = seq[j:j + args.wrap]
                    print(f"  {format_sequence(chunk, not args.no_seq_color)}")
            else:
                print(f"  {format_sequence(seq, not args.no_seq_color)}")
            
            # Plus line
            plus_content = plus[1:] if len(plus) > 1 else ''
            if plus_content:
                print(f"  {Colors.DIM}+{Colors.RESET}{format_header('+' + plus_content, '+')[len(Colors.DIM) + 1 + len(Colors.RESET):]}")
            else:
                print(f"  {Colors.DIM}+{Colors.RESET}")
            
            # Quality
            if args.wrap > 0:
                for j in range(0, len(qual), args.wrap):
                    chunk = qual[j:j + args.wrap]
                    qual_line = f"  {format_quality(chunk)}"
                    if args.raw_quality:
                        qual_line += f"  {Colors.DIM}{chunk}{Colors.RESET}"
                    print(qual_line)
            else:
                qual_line = f"  {format_quality(qual)}"
                if args.raw_quality:
                    qual_line += f"  {Colors.DIM}{qual}{Colors.RESET}"
                print(qual_line)
            
            print()  # Blank line between records
            
    except Exception as e:
        print(f"Error reading FASTQ file: {e}", file=sys.stderr)
        sys.exit(1)
    


if __name__ == '__main__':
    main()