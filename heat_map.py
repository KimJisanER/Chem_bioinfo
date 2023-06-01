import matplotlib.pyplot as plt
import pandas as pd
from tkinter import filedialog
import numpy as np
import seaborn as sns

# CSV 파일 선택
csv_file = filedialog.askopenfilename(title='Select the csv_file')
df = pd.DataFrame(pd.read_csv(csv_file, index_col=0))

# 그림 사이즈 지정
fig, ax = plt.subplots(figsize=(12, 8))

# 상삼각 행렬 마스크 생성
mask = np.zeros_like(df, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True


# mask = np.zeros_like(df, dtype=np.bool)
# mask[np.where(df == 0)] = True

heatmap=sns.heatmap(df,
            cmap = 'YlGnBu',
            annot=True,    # 실제 값을 표시한다
            mask=mask,     # 표시하지 않을 마스크 부분을 지정한다
            cbar_kws={"shrink": .5},  # 컬러바 크기 절반으로 줄이기
            vmin = -0.2, vmax=1)


# 플롯 다듬기
plt.title("5_HTR Distance Map")
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()

# 플롯 저장 및 출력
plt.savefig(csv_file + '_distancemap_nj_identity.png')
plt.show()

