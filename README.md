# fairnessclique


编译命令：g++ -std=c++11 main.cpp Graph.cpp RelativeFair.cpp baseline.cpp  -o main

## 弱公平团
./main 图.txt 属性.txt FairClique threshold order
### 不同order
基于节点排序： sorted
基于节点的度： degree
基于colorful degree： colorful
基于degeneracy: degeneracy


## 强公平团
./main 图.txt 属性.txt StrongClique threshold order 
### 不同order
基于节点排序： sorted
基于fairness degree： fairness
fairness degree的高维: heuristic


## 基于属性交替枚举的相对公平团
./main 图.txt 属性.txt RelativeStrong threshold order delta
### 不同order
基于节点排序： sorted
基于colorful degree： colorful
基于degeneracy: degeneracy
基于相对公平度： relative


## 基于弱公平团的相对公平团
./main 图.txt 属性.txt RelativeWeak threshold none delta