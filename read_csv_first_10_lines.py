import pandas as pd

def print_first_10_lines_of_csv(file_path):
    try:
        df = pd.read_csv(file_path)
        print(df.head(10))
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    print_first_10_lines_of_csv('output_method2_metrics.csv')
