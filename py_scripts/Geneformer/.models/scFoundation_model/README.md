---
language:
- en
---

# scFoundation

## Overview
scFoundation is a publicly available single-cell foundation model developed by researchers at BioMap, Tsingua University, and MBZUAI.

Additional details about model design and performance can be found in [Hao et al. 2024](https://www.nature.com/articles/s41592-024-02305-7).

The code can be found on the official [GitHub](https://github.com/biomap-research/scFoundation). 

## Usage Examples

This backbone can be adapted for arbitrary downstream tasks with ModelGenerator - please consult the [documentation](https://genbio-ai.github.io/ModelGenerator/) for more details. 

### Cell Type Classification
```
mgen fit --model SequenceClassification --model.backbone scfoundation --data CellClassificationDataModule --data.path <hf_or_local_path_to_your_dataset>
mgen test --model SequenceClassification --model.backbone scfoundation --data CellClassificationDataModule --data.path <hf_or_local_path_to_your_dataset>
```

## Limitations
* We only include scFoundation's encoder in ModelGenerator. The decoder is disabled. 

## Citation
```
@article{hao2024large,
  title={Large-scale foundation model on single-cell transcriptomics},
  author={Hao, Minsheng and Gong, Jing and Zeng, Xin and Liu, Chiming and Guo, Yucheng and Cheng, Xingyi and Wang, Taifeng and Ma, Jianzhu and Zhang, Xuegong and Song, Le},
  journal={Nature methods},
  volume={21},
  number={8},
  pages={1481--1491},
  year={2024},
  publisher={Nature Publishing Group US New York}
}
```
