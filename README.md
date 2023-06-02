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
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок6.jpg" width="400">
</p>
2. В 99% случаев интроны начинаются с последовательности GT (GU в РНК) и заканчиваются последовательностью AG. 
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок2.png" width="400">
</p>
3. В результате сплайсинга кодон может быть разделен интроном тремя способами. Следовательно, число нуклеотидов в экзонах или интронах НЕ ОБЯЗАТЕЛЬНО кратно трем.
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок5.png" width="400">
</p>
4. Гены, кодирующие белок, транскрибируются РНК-полимеразой II с образованием мРНК. Промоторы многих генов, транскрибируемых полимеразой II, содержат последовательность, подобную TATAA, на 25–30 нуклеотидов выше сайта начала транскрипции. 
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок3.png" width="400">
</p>
5. Консенсусная последовательность Козака (Консенсус Козака или последовательность Козака) представляет собой мотив нуклеиновой кислоты, который функционирует как сайт инициации трансляции белка в большинстве эукариотических транскриптов мРНК.

 *(gcc)gccRccATGG*

- Подчеркнутые нуклеотиды указывают кодон начала трансляции , кодирующий метионин, если рассматривать для ДНК, то ATG соответственно.
- Заглавные буквы обозначают высоко консервативные основания , т.е. последовательность «AUGG» постоянна или редко, если вообще меняется.  
- «R» означает, что в этом положении всегда наблюдается пурин (аденин или гуанин) (по Козаку чаще встречается аденин).
- Строчная буква обозначает наиболее распространенное основание в положении, где основание, тем не менее, может варьироваться.
- Последовательность в скобках (gcc) имеет неопределенное значение.
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок4.png" width="400">
</p>
6. На расстоянии приблизительно 300 п.н. от стоп кодона расположена консенсусная последовательность AATAAA.
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Рисунок7.jpg" width="400">
</p>

## Используемый датасет

Для определения начала и конца экзон-интронных последовательностей используем особые участки ДНК – нетранслируемые области (или нетранслируемые последовательности, НТП). Они не используются в качестве матрицы для синтеза белка и прилегают к транслируемой области (т.е. искомому участку) с обеих сторон: в начале – 5'НТП, в конце – 3'НТП. 

Для решения данной задачи создается 2 модели МО и, соответственно, 2 датасета: первая для определения начала экзон-интронной последовательности, вторая – для определения ее конца.
Используемый датасет включает в себя участки 5'НТП, 3'НТП и участки, не являющиеся ни 5'НТП, ни 3'НТП (межгенные области или транслируемые области). Причем для модели, определяющей начало искомого участка, класс "1" включает в себя 5'НТП, класс "0" – 3'НТП и участки, не являющиеся ни 5'НТП, ни 3'НТП. А для модели, определяющей конец искомого участка, наоборот, класс "1" включает в себя 3'НТП, а класс "0" – 5'НТП и участки, не являющиеся ни 5'НТП, ни 3'НТП.

Для снижения объема обрабатываемой информации датасет создается только из последовательности первой хромосомы человека.

Опишем процесс создания датасетов:
1. Выделение 5'НТП и 3'НТП последовательностей 
2. Выделение последовательностей, не являющихся ни 5'НТП, ни 3'НТП
3. Исключение пересечений последовательностей, не являющихся ни 5'НТП, ни 3'НТП, между собой и с 5' и 3'НТП последовательностями
4. Разделение последовательностей на классы "0" и "1" для двух моделей МО

## Оценка моделей 
<p class="text-justify">В качестве конечного результата мы хотим получать модель, которая может определять участки интереса в последовательности нуклеотидов ДНК.

Пример выходных данных: 
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/svm_eval_067.png" width="600">
</p>
<p class="text-justify"> На графике выше красный прямоугольник - область интереса, синяя кривая - предсказанная моделью вероятность.  
Если обозначить длину последовательности интереса за l, то модель смотрит на 2l нуклеотидов левее левой границы области и на 2l правее правой границы.  
Таким образом, при оценке модель рассматривает участок последовательности длинной 5l.   
</p>

<p class="text-justify"> Для каждого нуклеотида **x**, предсказанное значение **y** представляет собой вероятность наличия области интереса в окрестности **x**.  
Размер окрестности (окна) является гиперпараметром, в работе мы используем сразу несколько таких окон и суммируем предсказания модели в них.
</p>

<p class="text-justify">В идеале модель должна выдавать близкие к 1 значение в пределах красной области, и близкие к 0 значения за ее пределами.  
Ниже приведен пример с низкой точностью предсказания (т.к. в красной области значения предсказанной вероятности невысокие).   
</p>

<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/svm_eval_0.png" width="600">
</p>

Численно точность прогнозирования оценивается с помощью Dice метрики, которая обычно применяется для решения задач сегментации.
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/Dice.jpg" width="800">
</p>

В нашем случае:
- А - область интереса, площадь которой рассчитывается как 1.0 * l, где l - длина последовательности интереса (количество нуклеотидов)
- B - область предсказания, площадь которой представляет собой интеграл функции вероятности по всему рассматриваемому участку (5l)
- A∪B - область пересечения A и B, площадь которой представляет собой интеграл функции вероятности по области интереса (1l)

Для оценки модели мы используем **среднее значение и стандартное отклонение** распределения dice-метрик, рассчитанных на всех последовательностях класса 1.
<p align="center">
<img src="https://github.com/polilya/gene_predict_project/blob/main/figures/dice_dist.png" width="600">
</p>

## Решение задачи алгоритмами классического ML
В качестве модели для решения задачи поиска нетранслируемой последовательности был выбран SVM  
Всего обучалось **2 модели**:  
1-ая определяет последовательность **перед первой cds** в пределах первого экзона   
2-ая определяет последовательность **после последней cds** в пределах последнего экзона   

**Pipeline** состоит из 3 основных этапов:
- формирование датасета
- настройка гиперпараметров
- обучение и тестирование модели

Для подбора гиперпараметров SVM моделей используется отдельный скрипт **model_selection.py** в котором задействована функция RandomSearchCV из пакета scikiy-learn. 
После подбора, гиперпараметры автоматически вписываются в конфиг той модели, для которой осуществлялся подбор.

Экспериментально было выяснено, что выбор максимизируемой в процессе подбора гиперпараметров метрики **не оказывает влияния** на результирующие параметры распределения dice-метрик.


Зависимость распределения метрики от объема датасета (количества примеров класса 0)
|Size of zero class |Mean             |Std            |
|:----------------- |:---------------:|:-------------:|
| **1000**          | 0.33            | 0.06          |
| 2000              | 0.33            | 0.07          |
| 4000              | 0.32            | 0.09          |
| 8000              | 0.29            | 0.12          |
| 16000             | 0.27            | 0.14          |
| 32000             | 0.23            | 0.17          |

С увеличением объема датасета снижается среднее значение Dice-метрик.  
Вероятно, извлекаемые на этапе препроцессинга фичи не релевантны и большое количество примеров класса с похожими (на 1 класс) фичами путает модель.  


Зависимость распределения метрики от используемых окон (окрестностей нуклеотида)
|List of windows    |Mean             |Std            |
|:----------------- |:---------------:|:-------------:|
| 10                | 0.34            | **0.02**      |
| **25**            | **0.35**        | **0.03**      |
| 100               | **0.35**        | 0.10          |
| 200               | 0.28            | 0.14          |
| 500               | 0.22            | 0.16          |
| 3, 5, 7, 11       | 0.33            | **0.01**      |
| 10, 25            | 0.34            | **0.02**      |
| 10, 25, 50        | **0.35**        |" 0.03         |
| 50, 100, 200, 300 | 0.33            | 0.08          |
| 100, 200, 500     | 0.30            | 0.11          |
| 500, 1000, 1500   | 0.19            | 0.15          |

## Решение задачи с применением DL




## Сравнение результатов
