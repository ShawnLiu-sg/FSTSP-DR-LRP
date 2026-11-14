"""Convert one historical Solomon input to the C++ FSTSP_DR_LRP 1 format."""

import argparse
from decimal import Decimal
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--release-scale', required=True, type=Decimal,
                        help='scale multiplying the stored DELAY_TIME coefficients')
    parser.add_argument('--rendezvous-count', required=True, type=int,
                        help='number of nearest customer nodes eligible for replenishment')
    args = parser.parse_args()
    if not args.release_scale.is_finite() or args.release_scale <= 0:
        parser.error('--release-scale must be positive and finite')

    rows = [line.split() for line in args.source.read_text(encoding='utf-8-sig').splitlines()[1:] if line.strip()]
    depot, *customers = rows
    n = len(customers)
    if not 0 <= args.rendezvous_count <= n:
        parser.error('--rendezvous-count must be between zero and the number of customers')
    x, y = map(Decimal, depot[1:3])

    def distance(row):
        return ((Decimal(row[1]) - x)**2 + (Decimal(row[2]) - y)**2, int(row[0]))

    rendezvous = {int(row[0]) for row in sorted(customers, key=distance)[:args.rendezvous_count]}
    output = ['FSTSP_DR_LRP 1', f'{n} 0.1 0.1 4.5 11 1000', ' '.join(depot[1:3])]
    for row in customers:
        scale_factor = Decimal(row[6]) if len(row) == 7 else Decimal(0)
        node = int(row[0])
        output.append(' '.join([*row[:6], str(scale_factor * args.release_scale),
                                str(int(scale_factor > 0)), str(int(node in rendezvous))]))
    args.destination.write_text('\n'.join(output) + '\n', encoding='ascii')


if __name__ == '__main__':
    main()
