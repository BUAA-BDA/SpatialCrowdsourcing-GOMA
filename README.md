GOMA: Two-sided Online Micro-Task Assignment in Spatial Crowdsourcing
========================================================================

This repository stores the source code of the solutions to the problem called GOMA in the following papers.
If you use our source code or dataset, please consider citing our papers.

[1] **Online mobile Micro-Task Allocation in spatial crowdsourcing.**
*Yongxin Tong, Jieying She, Bolin Ding, Libin Wang, Lei Chen.* ICDE 2016: 49-60. [link](https://doi.org/10.1109/ICDE.2016.7498228)
 
[2] **Two-sided Online Micro-Task Assignment in Spatial Crowdsourcing.**
*Yongxin Tong, Yuxiang Zeng, Bolin Ding, Libin Wang, Lei Chen.* IEEE Transactions on Knowledge and Data Engineering, 2019. [link](https://doi.org/10.1109/TKDE.2019.2948863)
 


Usage of the algorithms
---------------

### Environment

gcc/g++ version: 7.4.0 

OS: Ubuntu

### Compile the algorithms

cd algorithm && make all


### Run the algorithms

./TGOA-OP ./realData/EverySender\_cap1/data\_00.txt

./Greedy ./realData/EverySender\_cap1/data\_00.txt

./TGOA ./realData/EverySender\_cap1/data\_00.txt

./TGOA-Greedy ./realData/EverySender\_cap1/data\_00.txt

./Ext-GRT ./realData/EverySender\_cap1/data\_00.txt

./OPT ./realData/EverySender\_cap1/data\_00.txt

Description of the datasets
---------------

### Environment

Python: 2.7

### Synthetic dataset

dataset/synthetic: a sample of our synthetic dataset (\#2)

dataset/genDataSynthetic.py: a script to generate the synthetic datasets

Please refer to genDataSynthetic.py for the format of the dataset.

### Real dataset

dataset/real/EverySender*: includes the datasets of EverySende 

dataset/real/gMission*: incldues the datasets of gMission

Please refer to the source code for the format of the dataset.



Contact
------------
- Yuxiang Zeng: yzengal@cse.ust.hk
- Yongxin Tong: yxtong@buaa.edu.cn
