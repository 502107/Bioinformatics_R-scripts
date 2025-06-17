import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as file:
    lines = file.readlines()
    start_line = next(i for i, line in enumerate(lines) if line.startswith('GO'))
    data = pd.read_csv(input_file, delimiter='\t', skiprows=start_line)

data = data[['GO','name', 'p_fdr_bh', 'NS', 'study_count', 'ratio_in_study']].rename(columns={
    'name': 'term',
    'p_fdr_bh': 'adjusted_p_value',
    'NS': 'source',
    'study_count': 'intersection_size'
})

data = data[data['adjusted_p_value'] < 1]

# def count_leading_dots(term):
#     return len(term) - len(term.lstrip('.'))
# data['dot_count'] = data['GO'].apply(count_leading_dots)
# print(data['dot_count'])
# print(data['GO'])
# data = data[data['dot_count'] > 0]
# data = data.drop(columns=['dot_count'])

data.to_csv(output_file, index=False)