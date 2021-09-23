# Cancer Somatic Variant Classification Tool (CSVC)

依据肿瘤患者个人信息和变异数据，进行变异分类和靶向药物注释。

# Installation

无需进行安装

# Usage

使用CSVC.py来运行CSVC，进行靶向药物推荐,代码如下:

```
python3 CSVC.py -i <input_somatic> -n <sample_name> -o <outputdir>
```

CSVC_FUSCC只接受三个参数：

1.Somatic变异的输入文件：需要给出相对或绝对路径，e.g. ./test/1372.vcf

2.样本名	e.g. 1372

3..输出目录 ：需要给出相对或绝对路径，并且该路径需要存在 e.g.  ./ouput



# Version

CSVC_20210922