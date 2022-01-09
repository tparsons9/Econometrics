# Econometrics
Collection of Econometrics texts and syllabi, along with attempted problem sets. 
The following is a curated reading list taken from various courses syllabi (GWU PPPA 6022, MIT 14.32 (2007), MIT 14.36(2013), and BU EC508-A1 (2021)). 
Problem sets and lecture notes are borrowed from MIT 14.382 (2017)

-- **Note:** This is my personal attempt folder for learning econometric methods. All credit goes to Prof(s): Angrist, Hausman, Smith, Forneron, and Chernozhukov for their contributions to course materials. This is merely a collection of materials meant to better my own learning. This repository is public to allow other interested learners an outline combining various texts and levels of difficulty. My general understanding may change, but I have included texts from Wooldridge (Undergraduate - Foundations of Econometrics), Angrist & Pischke (Graduate - Advanced Econometrics), and Hayashi (Graduate/PhD - Theoretical Supplement). I have also included selected papers taken from the referenced syllabi to serve as supplements to the reading material. 

As my research aims take a marketing focus, I have also included selected readings from Prof. Avi Goldfarb and his PhD Course RSM 3052 at the University of Toronto. I am deeply grateful for his generosity in providing this reading list. Marketing-related papers will be indicated by a "**+**", and will be related to the Econometric course material at hand. 

Finally, I will include my interpretations of the assignments with the provided data from Angrist/Chernozhukov. Since assignment answers are not provided, I make no claim to the veracity of the results. 

As I iterate through this "course," I will find things that may not fit, or doesn't match the flow of learning that I am trying to achieve. This is afterall a mixture of individual ideas on the best approach to learning econometric methods, and I'm sure I'll find ways along the way to make the course more cohesive and complete. I ask that you bear with me (is there anybody there?) as I work through this process. As for where I'm turning my attention after completion, I'll be auditing Grant McDermott's fantastic course [Environmental Economics and Data Science (EC607)](https://github.com/uo-ec607), hosted on GitHub. Thanks, and happy learning!

All references are included at the end of this ReadMe file. (Papers are not properly cited in Outline-- they were copied from various sources). PDF's will be in their respective topic folders, along with problem sets and data. My attempt at the problem sets will be in a separate branch. 

## Course Outline 

### Statistics Review 
* **Text** 
  - Wooldridge -- Apendix B & C
* **Lecture Notes** 
  - None
* **Papers** 
  - Zsohar, Peter. (2012). Short Introduction to the Generalized Method of Moments. Statisztikai szemle: a Magyar Központi      Statisztikai Hivatal folyóirata. 90. 150-170. 
  - Burtless, Gary. "Are Targeted Wage Subsidies Harmful? Evidence from a Wage Voucher Experiment." Industrial and Labor Relations Review 39 (October 1985): 105-111.
* **Assignments** 
  -  14.32 -- Review Problem Set & Problem Set 1  
  
### Simple Linear Regression
* **Text**
  - Wooldridge: 1-3
  - Angrist & Pischke: 1, 2, 3.1-3.2
  - Hayashi: 1-3
* **Lecture Notes**
  - Chernozhukov: L1. Least Squares, Adaptive Partialling-Out, Simultaneous Interference
* **Papers** 
  - None
* **Assignments** 
  -  Homework 1 (All assignments are from 2017 MITOCW 14.32, unless stated otherwise) 
  -  14.32 -- Problem Set II

### Multivariate Regression, Variability, and Instrumental Variables
* **Text** 
  - Wooldridge: 4-7, 15 & 19 
  - Angrist: Ch. 4
  - Hayashi: Ch. 4
* **Lecture Notes** 
  - Chernozhukov: L2. Structural Equation Models and IV
  - Chernozhukov: L3. Structural Equation Models and GMM
* **Papers** 
  - Hansen, Hausman and Newey (2008): “Estimation with Many Instrumental Variables,”
Journal of Business and Economic Statistics.
  -Hahn and Hausman (2003), “IV Estimation with Valid and Invalid Instruments,”
mimeo
  - Moorthy, K.S. 1993. "Theoretical Modeling in Marketing." Journal of Marketing, 57 (April): 92- 106. **+**
  - Shmueli, G. 2010. “To Explain or To Predict?” Statistical Science, 25 (3): 289-310. **+**
* **Assignments** 
  -  14.32 -- Problem Set III-IV 
  -  Homework 2 and 3


### Nonlinear Regression 
* **Text** 
  - Angrist: 3.3-3.4 & 7
  - Hayashi: 7, 8.3
  - [Non-linear Least Squares](https://towardsdatascience.com/a-guide-to-building-nonlinear-least-squares-nls-regression-models-310b97a7baeb)
* **Lecture Notes** 
  - Chernozhukov: L4. Euler Equations, Nonlinear GMM, and Other Adventures
  - Chernozhukov: L5. Bootstrapping
  - Chernozhukov: L6. Nonlinear and Binary Regression, Predictive Effects, and M-Estimation
* **Papers** 
  - Chay, K., and J. Powell. 2001. "Semiparametric Censored Regression Models." Journal of Economic Perspectives, 15: 29-42.
* **Assignments** 
  -  Homeworks 4 & 5

### Panel Data 
* **Text** 
  - Wooldridge: 10 
  - Hayashi: 5 
  - Angrist: 5
* **Lecture Notes** 
  - Chernozhukov: L8. Linear Panel Data Models Under Strict and Weak Exogeneity 
  - Chernozhukov: L10. Nonlinear Panel Data 
* **Papers** 
  - Hausman, J., and W. Taylor. 1981. "Panel Data and Unobservable Individual Effects." Econometrica, 49(6): 1377-1398. 
  - Griliches, Z., and J. Hausman. 1986. "Errors in Variables in Panel Data." Journal of Econometrics, 31: 93-118. 
* **Assignments** 
  -  None


### Heteroskedasticity, Time Series, and Serial Correlation
* **Text** 
  - Wooldridge: 8, 11.3, 11.4, 12, 13
  - Hayashi: 6, 9 (10?) 
* **Lecture Notes** 
  - [Heteroskedasticity and Robust Estimators](http://www3.grips.ac.jp/~yamanota/Lecture_Note_9_Heteroskedasticity)
  - [Serial Correlation](https://www3.nd.edu/~rwilliam/stats2/l26.pdf)
* **Papers** 
  - Graddy, K. 1995. "Testing for Imperfect Competition at the Fulton Fish Market." RAND Journal of Economics 26, no. 1: 75-92.
* **Assignments** 
  -  14.32 Problem Set V
 
### Simultaneous Equations and Applied Models
* **Text** 
  - Wooldridge: 16
* **Lecture Notes** 
  - Chernozhukov: L7. Distribution Regression and Counterfactual Analysis
* **Papers** 
  - Angrist, J., G. Imbens, and K. Graddy. "The Interpretation of Instrumental Variables Estimators in Simultaneous Equations Models with an Application to the Demand for Fish." Review of Economic Studies 67, no. 3 (July 2000): 499-527.
  - Hausman, Hall, and Griliches (1984): “Econometric Models for Count Data with an
Application to the Patents-R & D Relationship,” Econometrica. 
  - Han and Hausman (1990): “Flexible Parametric Estimation of Duration and Competing
Risk Models,” Journal of Applied Econometrics
  - McFadden and Train (2000): “Mixed MNL Models for Discrete Response,” Journal of
Applied Econometrics
  - [Train (2003): “Discrete Choice Methods with Simulation”](https://eml.berkeley.edu/books/choice2.html)
  - Burda, Harding and Hausman (2008): “A Bayesian Mixed Logit-Probit Model for Multinomial Choice,” Journal of Econometrics.
* **Assignments** 
  -  Homework 6
  -  14.32 Problem Set VI

## Other Useful Sources: 
* [EC508 Notes & Review](https://alex-hoagland.github.io/teaching/2019-ec508)
* Chernozhukov L9. GMM Under Moderately High Dimensions
* Chernozhukov L11. Inference for High-Dimensional Sparse Econometric Models
* Chernozhukov L12. Treatment Effects

## References 
### Texts
* Wooldridge, Jeffrey M., 1960-. Introductory Econometrics : a Modern Approach. Mason, Ohio :South-Western Cengage Learning, 2006.
* Angrist, J. D., and Jorn-Steffen Pischke. 2008. Mostly Harmless Econometrics. Princeton, NJ: Princeton University Press.
* Hayashi, Fumio. Econometrics. Princeton :Princeton University Press, 2000.

### Syllabi 
* Smith, Jonathan. "Econometrics for Policy Research II," (Online, George Washington University, 2012).
* Angrist, Joshua. "14.32 Econometrics," (Online, Massachussets Institute of Technology, 2007). 
* Hausman, Jerry. "14.36 Undergraduate Advanced Econometrics," (Online, Massachussets Institute of Technology, 2013).
* Forneron, Jean-Jacques. "EC508-A1: Econometrics," (Online, Boston University, 2021).
* Goldfarb, Avi. "RSM 3052 Course Outline," (Online, University of Toronto, 2022). 

### Papers

- Chay, Kenneth Y., and James L. Powell. “Semiparametric Censored Regression Models.” Journal of Economic Perspectives 15, no. 4 (December 2001): 29–42. https://doi.org/10.1257/jep.15.4.29.
- Graddy, Kathryn. “Testing for Imperfect Competition at the Fulton Fish Market.” The RAND Journal of Economics 26, no. 1 (1995): 75. https://doi.org/10.2307/2556036.
- Griliches, Zvi, and Jerry A. Hausman. “Errors in Variables in Panel Data.” Journal of Econometrics 31, no. 1 (February 1, 1986): 93–118. https://doi.org/10.1016/0304-4076(86)90058-8.
- Han, Aaron, and Jerry A. Hausman. “Flexible Parametric Estimation of Duration and Competing Risk Models.” Journal of Applied Econometrics 5, no. 1 (1990): 1–28.
- Hansen, Christian, Jerry Hausman, and Whitney Newey. “Estimation With Many Instrumental Variables.” Journal of Business & Economic Statistics 26, no. 4 (October 1, 2008): 398–422. https://doi.org/10.1198/073500108000000024.
- Hausman, Jerry A., and William E. Taylor. “Panel Data and Unobservable Individual Effects.” Econometrica 49, no. 6 (November 1981): 1377. https://doi.org/10.2307/1911406.
- Hausman, Jerry, Bronwyn H. Hall, and Zvi Griliches. “Econometric Models for Count Data with an Application to the Patents-R & D Relationship.” Econometrica 52, no. 4 (1984): 909–38. https://doi.org/10.2307/1911191.
- Imbens, Guido, Joshua Angrist, and Kathryn Graddy. “The Interpretation of Instrumental Variables Estimators in Simultaneous Equations Models with an Application to the Demand for Fish.” Review of Economic Studies 67, July (2000): 499–527.
- Moorthy, K. Sridhar. “Theoretical Modeling in Marketing.” Journal of Marketing 57, no. 2 (April 1993): 92. https://doi.org/10.2307/1252029.


