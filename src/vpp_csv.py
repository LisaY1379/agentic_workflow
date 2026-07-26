import csv
import sys
from math import gcd, isqrt


# --- Pomerance triple verification function ---
def pp_verify(p, A, x0):
    if p < 5 or p % 2 == 0:
        return False

    q = isqrt(p)  # p not prime ==> prime divisor < q
    k = (q + 1 + isqrt(4 * q)).bit_length()  # least k such that 2^k > q + 1 + 2*sqrt(q)

    if gcd(A * A - 4, p) != 1:
        return False  # singular

    X, Z = x0 % p, 1
    Zprev = None

    # precompute C = (A + 2)/4 mod p for doubling formula
    C = ((A + 2) * ((p + 1) // 4 if p % 4 == 3 else (3 * p + 1) // 4)) % p
    for i in range(k):
        Zprev = Z
        U, V = (X + Z) * (X + Z) % p, (X - Z) * (X - Z) % p
        W = U - V
        X, Z = U * V % p, W * (V + C * W) % p

    # if p|Z and (Zprev,p)=1 then (p,A,x0) is valid
    return Z % p == 0 and gcd(Zprev, p) == 1


# --- CSV batch processing and verification logic ---
def verify_csv(file_path):
    print(f"Reading and verifying file: {file_path}")

    total_count = 0
    skipped_count = 0
    false_count = 0
    false_rows = []

    with open(file_path, mode='r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)

        # Check if required header columns exist
        required_cols = {'prime', 'A', 'x0'}
        if not required_cols.issubset(set(reader.fieldnames)):
            print(f"Error: CSV must contain the following header columns: {required_cols}")
            return

        for row_idx, row in enumerate(reader, start=2):  # Start at row 2 (row 1 is header)
            total_count += 1

            # --- Skip invalid or broken data rows ---
            try:
                # Check for empty fields
                if not row['prime'] or not row['A'] or not row['x0']:
                    raise ValueError("Contains empty values")

                p = int(row['prime'].strip())
                A = int(row['A'].strip())
                x0 = int(row['x0'].strip())

                # Run verification
                is_valid = pp_verify(p, A, x0)

                if not is_valid:
                    false_count += 1
                    false_rows.append((row_idx, p, A, x0))

            except (ValueError, TypeError, KeyError) as e:
                # Catch parsing errors, log them, and skip the row
                skipped_count += 1
                print(f"[SKIP/FAILED] Line {row_idx} failed to parse ({e}), skipped. Data: {row}")
                continue

    # --- Summary output ---
    valid_tested = total_count - skipped_count
    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY:")
    print(f"Total rows read: {total_count}")
    print(f"Failed to parse (Skipped): {skipped_count}")
    print(f"Valid rows tested: {valid_tested}")
    print(f"Passed (True): {valid_tested - false_count}")
    print(f"Failed (False): {false_count}")

    if false_count > 0:
        print("\nResult: Found data that evaluated to False! ❌")
        print("\nFirst 10 failed rows (Line Number, prime, A, x0):")
        for row in false_rows[:10]:
            print(f"  Line {row[0]}: p={row[1]}, A={row[2]}, x0={row[3]}")
        if len(false_rows) > 10:
            print(f"  ... and {len(false_rows) - 10} more row(s)")
    else:
        print("\nResult: All valid rows passed verification (no False detected)! ✅")
    print("=" * 40)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <your_csv_file.csv>")
    else:
        verify_csv(sys.argv[1])