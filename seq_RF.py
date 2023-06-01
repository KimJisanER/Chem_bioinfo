import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tkinter import filedialog
import os
from BLOSUM62 import embed_sequence

# Prompt the user to select the fasta file
fasta_file = filedialog.askopenfilename(title='Select the fasta file')
target = input('input_accsesion_to_work_with:',)
pdb_file = filedialog.askopenfilename(title='Select the pdbfile')

# Open the file
with open(fasta_file, 'r') as f:
    lines = f.readlines()

seqs = []
headers = []
header_positions = []  # Store the positions of headers satisfying the condition
for i, line in enumerate(lines):
    if line.startswith('>'):
        headers.append(line.strip().upper())
        seqs.append('')
        if '>SP' in line and 'HUMAN' in line and 'HTR' in line:
            header_positions.append(i)  # Store the position of the header
    else:
        seqs[-1] += line.strip()

# Define the dataset dictionary
dataset = []

for header, seq in zip(headers, seqs):
    sequence_data = {
        'sequence': list(seq),
        'class': 0  # Initialize as non-5-HTR receptor
    }

    if '5-hydroxytryptamine receptor'.upper() in header or 'HTR' in header:
        sequence_data['class'] = 1

        # if '5-hydroxytryptamine receptor 5A'.upper() in header or 'HTR5A' in header\
    #         or '5-hydroxytryptamine receptor 7'.upper() in header or 'HTR7' in header \
    #         or '5-hydroxytryptamine receptor 1A'.upper() in header or 'HTR1A' in header \
    #         or '5-hydroxytryptamine receptor 1B'.upper() in header or 'HTR1B' in header \
    #         or '5-hydroxytryptamine receptor 1D'.upper() in header or 'HTR1D' in header:
    #     sequence_data['class'] = 1  # Set as 5-HTR receptor

    # if '5-hydroxytryptamine receptor 5A'.upper() in header or 'HTR5A' in header \
    #         or '5-hydroxytryptamine receptor 2B'.upper() in header or 'HTR2B' in header \
    #         or '5-hydroxytryptamine receptor 2C'.upper() in header or 'HTR2C' in header \
    #         or '5-hydroxytryptamine receptor 1A'.upper() in header or 'HTR1A' in header \
    #         or '5-hydroxytryptamine receptor 1B'.upper() in header or 'HTR1B' in header \
    #         or '5-hydroxytryptamine receptor 1D'.upper() in header or 'HTR1D' in header:
    #     sequence_data['class'] = 1  # Set as 5-HTR receptor
    #
    # if '5-hydroxytryptamine receptor 2A'.upper() in header or 'HTR2A' in header \
    #         or '5-hydroxytryptamine receptor 7'.upper() in header or 'HTR7' in header:
    #     sequence_data['class'] = 2  # Set as 5-HTR receptor


    # if '5-hydroxytryptamine receptor 6'.upper() in header or 'HTR6' in header\
    #         or '5-hydroxytryptamine receptor 4'.upper() in header or 'HTR4' in header \
    #         or '5-hydroxytryptamine receptor 1B'.upper() in header or 'HTR1B' in header \
    #         or '5-hydroxytryptamine receptor 1D'.upper() in header or 'HTR1D' in header:
    #     sequence_data['class'] = 1  # Set as 5-HTR receptor
    #
    # if '5-hydroxytryptamine receptor 2A'.upper() in header or 'HTR2A' in header \
    #         or '5-hydroxytryptamine receptor 2B'.upper() in header or 'HTR2B' in header \
    #         or '5-hydroxytryptamine receptor 2C'.upper() in header or 'HTR2C' in header \
    #         or '5-hydroxytryptamine receptor 1A'.upper() in header or 'HTR1A' in header \
    #         or '5-hydroxytryptamine receptor 7'.upper() in header or 'HTR7' in header:
    #     sequence_data['class'] = 2  # Set as 5-HTR receptor
    #     print(header)

    dataset.append(sequence_data)

    if target.upper() in header:
        target_seq = seq

filtered_dataset = [data for data in dataset if len(data['sequence']) == len(dataset[0]['sequence'])]
print(len(filtered_dataset))

# bindsite_list=[]
# for i in range(2035,2040):
#     bindsite_list.append(i)
# for i in range(2056,2060):
#     bindsite_list.append(i)
# for i in range(2319,2324):
#     bindsite_list.append(i)
# for i in range(2328,2338):
#     bindsite_list.append(i)
# for i in range(2536,2595):
#     bindsite_list.append(i)
# for i in range(2608,2957):
#     bindsite_list.append(i)
# for i in range(5105,5112):
#     bindsite_list.append(i)
# for i in range(5128, 5133):
#     bindsite_list.append(i)
# for i in range(5301,5350):
#     bindsite_list.append(i)


bindsite_list = [i for i in range(18940)]

for data in filtered_dataset:
    updated_sequence = ''
    for i in range(len(data['sequence'])):
        if i in bindsite_list:
            updated_sequence += data['sequence'][i]
    data['sequence'] = list(updated_sequence)
    print(len(updated_sequence))

sequences = [data['sequence'] for data in filtered_dataset]
classes = [data['class'] for data in filtered_dataset]

# sequence_encoded = []
# for sequence in sequences:
#     encoded_sequence = [ord(aa) - ord('A') for aa in sequence]
#     sequence_encoded.append(encoded_sequence)
#
# print(sequence_encoded[0])

sequence_encoded = []
for sequence in sequences:
    encoded_sequence = embed_sequence(sequence)
    sequence_encoded.append(encoded_sequence)

sequence_encoded = np.array(sequence_encoded)
sequence_encoded = sequence_encoded.reshape(sequence_encoded.shape[0], -1)
print(sequence_encoded[0])
print(len(sequences[0]))
print(len(sequence_encoded[0]))

num_iterations = 100

# 반복적으로 학습 및 예측 수행
for i in range(num_iterations):

    print(i)
    X_train, X_test, y_train, y_test = train_test_split(sequence_encoded , classes, test_size=0.1, stratify=classes)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1, stratify=y_train)

    # 랜덤 포레스트 분류 모델 생성
    # rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_model = RandomForestClassifier(n_estimators=100)

    # feature importance를 저장할 리스트
    importance_scores = np.zeros(len(X_train[0]))

    # 랜덤 포레스트 모델 학습
    rf_model.fit(X_train, y_train)

    y_val_pred = rf_model.predict(X_val)
    val_accuracy = accuracy_score(y_val, y_val_pred)
    print("Validation Accuracy:", val_accuracy)

    y_test_pred = rf_model.predict(X_test)
    test_accuracy = accuracy_score(y_test, y_test_pred)
    print("Test Accuracy:", test_accuracy)

    # feature importance 계산
    importances = rf_model.feature_importances_
    # importance_scores에 추가
    importance_scores += importances

target_dict = {}
aa_count = 1
bindsite_count = 0
for i in range(len(target_seq)):
    if target_seq[i] == '-' and i in bindsite_list:
        bindsite_count += 1
        pass
    elif target_seq[i] != '-' and i in bindsite_list:
        target_dict[bindsite_list[bindsite_count]] = aa_count
        aa_count += 1
        bindsite_count += 1
    elif target_seq[i] == '-' and i not in bindsite_list:
        pass
    else:
        aa_count += 1
        # target_dict[i] = importance_scores[i]*10000

print(target_dict)

aa_count = 0
coloring_dict = {}
non_coloring_dict ={}
for i in range(len(bindsite_list)):
    if target_seq[bindsite_list[i]] == '-':
        if bindsite_list[i] not in target_dict.keys():
            non_coloring_dict[bindsite_list[i]+1] = importance_scores[i] * 1000
    else:
        if bindsite_list[i] in target_dict.keys():
            coloring_dict[target_dict[bindsite_list[i]]] = importance_scores[i] * 1000

print(non_coloring_dict)
print(coloring_dict)


aa_count = 0
coloring_dict = {}
non_coloring_dict ={}
for i in range(len(bindsite_list)):
    if target_seq[bindsite_list[i]] == '-':
        if bindsite_list[i] not in target_dict.keys():
            importance_score_24 = 0
            for j in range(24):
                importance_score_24 += importance_scores[24*i + j]
            non_coloring_dict[bindsite_list[i]+1] = importance_score_24 * 1000
    else:
        if bindsite_list[i] in target_dict.keys():
            importance_score_24 = 0
            for j in range(24):
                importance_score_24 += importance_scores[24 * i + j]
            coloring_dict[target_dict[bindsite_list[i]]] = importance_score_24 * 1000

print(non_coloring_dict)
print(coloring_dict)


def get_importance_pdb(pdb_file, coloring_dict):

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    for i in range(len(lines)):
        if lines[i].startswith('ATOM'):
            resi_num = int(lines[i][22:26])
            b_factor = float(lines[i][60:66])
            importance = 0
            if resi_num in coloring_dict.keys():
                importance = coloring_dict[resi_num]
            else: pass

            lines[i] = lines[i][:60] + f'{importance:6.2f}' + lines[i][66:]

        elif lines[i].startswith('HETATM'):
            importance = 0

            lines[i] = lines[i][:60] + f'{importance:6.2f}' + lines[i][66:]

    file_path = os.path.join(f'{os.path.basename(pdb_file)[:-4]}_HTR_all_530.pdb')
    with open(file_path, 'w') as f:
        f.writelines(lines)

get_importance_pdb(pdb_file, coloring_dict)
