# Разработка пайплайна проекта "Предсказание генов"
Предсказание генов – определение кодирующих и регуляторных последовательностей ДНК в геноме. Кодирующие последовательности - последовательности, несущие информацию о строении белков. Регуляторные последовательности – вспомогательные последовательности в цепочке преобразования последовательности ДНК в белок.

Цель проекта: обучение моделей ML/DL для определения экзон-интронных последовательностей.

## Обзор
У большинства организмов гены состоят из ДНК, где определенная последовательность ДНК определяет функцию гена. Ген транскрибируется с ДНК на РНК, которая может либо быть некодирующей (ncRNA) или промежуточным звеном (mRNA), траслируемым впоследствии в белок. Каждая из этих стадий контролируется определенным элементами последовательности или области внутри гена. Каждый функциональный  ген, следовательно, требует нескольких обязательных элементов последовательности. Это включает в себя последовательность, которая фактически кодирует функциональный белок, а также несколько областей регуляторной последовательности.

На рисунке ниже представлена структура кодирующего белок гена эукариотов. 

- The open read frame (ORF) представляет собой область, состоящую из экзонов и интронов, которая траскрибируется сначала на pre-mRNA, а затем в процессе слайсинга происходит "вырезание" интронных областей и образование конечной, готовой к трансляции mRNA. 

- Регуляторные последовательсноти расположены по обеим сторонам от ORF. 5'UTR и 3'UTR - это нетранслируемые участки, содержащие регуляторные элементы.
<p align="center">
   <img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок1.png" width="400">
</p>
### *Закономерности в строении последовательностей*

1. ATG - стартовый кодон, с которого начинается процесс транскрипции, TAA, TAG и TGA - стоп кодон, прерывающие данный процесс.

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок6.jpg" width="400">

2. В 99% случаев интроны начинаются с последовательности GT (GU в РНК) и заканчиваются последовательностью AG. 

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок2.png" width="400">

3. В результате сплайсинга кодон может быть разделен интроном тремя способами. Следовательно, число нуклеотидов в экзонах или интронах НЕ ОБЯЗАТЕЛЬНО кратно трем.

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок5.png" width="400">

4. Гены, кодирующие белок, транскрибируются РНК-полимеразой II с образованием мРНК. Промоторы многих генов, транскрибируемых полимеразой II, содержат последовательность, подобную TATAA, на 25–30 нуклеотидов выше сайта начала транскрипции. 

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок3.png" width="400">

5. Консенсусная последовательность Козака (Консенсус Козака или последовательность Козака) представляет собой мотив нуклеиновой кислоты, который функционирует как сайт инициации трансляции белка в большинстве эукариотических транскриптов мРНК.

    *(gcc)gccRccATGG*

- Подчеркнутые нуклеотиды указывают кодон начала трансляции , кодирующий метионин, если рассматривать для ДНК, то ATG соответственно.
- Заглавные буквы обозначают высоко консервативные основания , т.е. последовательность «AUGG» постоянна или редко, если вообще меняется.  
- «R» означает, что в этом положении всегда наблюдается пурин (аденин или гуанин) (по Козаку чаще встречается аденин).
- Строчная буква обозначает наиболее распространенное основание в положении, где основание, тем не менее, может варьироваться.
- Последовательность в скобках (gcc) имеет неопределенное значение.

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок4.png" width="400">

6. На расстоянии приблизительно 300 п.н. от стоп кодона расположена консенсусная последовательность AATAAA.

<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок7.jpg" width="400">

## Используемый датасет

## Решение задачи алгоритмами классического ML

## Решение задачи с применением DL

## Сравнение результатов
