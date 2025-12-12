## Relatório - Terceiro trabalho de Álgebra Linear
##### Escrito por Daniel Birenbaum, 11/12/2025
##### github.com/danielbirenbaum

### Definindo a função $dij$ e provando o rank de $G$
Dada a função $d_{ij} = f(\cdot)$, que representa os valores nas posições $i$ e $j$ da matriz $D$ de distâncias

$d_{ij} = f(x_i,x_j) = \sqrt{(x_{1i} - x_{1j})^2 + (x_{2i} - x_{2j})^2}$

Define-se $D_f$, $C_f$ e $I_f$ como domínio, contradomínio e conjunto imagem da função $f$.

$D_f = \{x_{1i},x_{1j},x_{2i},x_{2j} \in \mathbb{R};\ x = (x_1,x_2)\in \mathbb{R^2} \}$

$C_F =\mathbb{R_{\ge0}}$

$I_f = \{d_{ij} = f(x_i,x_j);\ d_{ij} \in [0,\infty] \}$

Se $x_{i} = [x_{1i} \ \ x_{2i}]^T$ são os elmenetos da matriz $X = [x_1 \ \ x_2 \ \ x_3 \ \ ... \ \ x_{27}]$, $X$ é $2\times27$, $X^T$ é $27\times2$. Logo $X^TX$ é quadrada de dimensão $n\times n$ com $n = 27$.
Define-se $G = X^T X$, perceba:

```math

G= X^TX \Rightarrow G^T = (X^TX)^T \Rightarrow G^T = X^TX = G

```



Entende-se que $G = G^T$, portanto simétrico.

Percebe-se também, que a matriz $X = [x_1 \ \ x_2 \ \ x_3 \ \ ... \ \ x_{27}]$ possui _rank_ $r(X)= 2$ pois $x_i = [x_{1i} \ \ x_{2i}]^T$. O span da matriz $X$, portanto, é no máximo, o $\mathbb{R}^2$, visto que se houver duas colunas linearmente independentes, será possível alcançar o $\mathbb{R}^2$. 

Sabe-se que o $r(X)= r(X^TX)$:
```math
Xv = 0 \Rightarrow X^TXv = 0\Rightarrow (X^TXv)^T = 0 \Rightarrow v^TX^TX = 0 \Rightarrow v^TX^TXv = 0 
```
```math
\Rightarrow (Xv)^T(Xv) = 0 \Rightarrow \parallel Xv \parallel^2 = 0
```

$X^TX$ possui o mesmo espaço nulo que X. 

```math
ker(X) = ker(X^TX) 
```
```math
r(X) + ker(X) = 27 
```
```math
r(X^TX) + ker(X^TX) = 27
```
```math
r(X) = r(X^TX)
```

Portanto, $G$ possui _rank_ igual ao da matriz $X$. $G$ é simétrica, e possui _rank_ de no máximo 2, logo terá no máximo 2 autovalores não nulos que formem, com os autovetores, o espaço coluna.

### Análise do operador linear $P$

Sabe-se que $P$ é idempotente, $P = P^2$, seja $\lambda$ e $v$ um autovalor e um autovetor de $P$:
```math
Pv = \lambda v \Rightarrow P^2v = P \lambda v \Rightarrow P^2v = \lambda^2v 
```
```math
\lambda v = \lambda^2v 
```
```math
\lambda^2 - \lambda = 0 
```
```math
\lambda_1 = 0 \ \ \ \ \lambda_2 = 1
```

Logo, somente é possível que os autovalores de $P$ assumam $0$ ou $1$. Para definir o autovetor a $\lambda = 0$, usa-se a definição de $P$, seja $u = [1 \ \ 1 \ \ 1 \ ...\  1]^T$:
```math
Pu = u - \frac{\sum_{l=1}^n u_l}{n} \cdot \mathbf{1} 
```
```math
\frac{\sum_{l=1}^n u_l}{n} = 1
```
```math
Pu = u - \mathbf{1}
```
```math
Pu = 0
```

Portanto, o $u = [1 \ \ 1 \ \ 1 \ ...\  1]^T$ é autovetor referente a $\lambda = 0$, o que é notório se considerarmos que a matriz $P$ projeta vetores $u$ ao espaço ortogonal ao espaço gerado pelo vetor de **1**'s.

Os vetores $v$ cujos produtos internos com **1** sejam zero, devem estar no espaço coluna de $P$, portanto:
```math
\mathbf{1}^Tv = 0 
```
```math
Pv = (I - \frac{\mathbf{1}\mathbf{1}^T}{n})v  
```
```math
Pv = v - \frac{\mathbf{1}(\mathbf{1}^Tv)}{n} 
```
```math
Pv = v
```

Logo vetores do tipo $v_{i=1} = [-1 \ \ 1 \ \ 0 \ \ ... \ \ 0]^T$, $v_{i=2} = [-1 \ \ 0 \ \ 1 \ \ ... \ \ 0]^T$,etc... Serão autovetores de $P$ para $\lambda = 1$. Os resultados são coerentes com aqueles computados no arquivo .py.

Seja $A = \frac{\mathbf{1}\mathbf{1}^T}{n}$ a matriz de projeção que induz $P$:
```math
A = VLV^{-1} 
```
```math
L =
\begin{bmatrix}
0 & 0 & \cdots &0 \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{bmatrix}
\\[1em]
```
```math
V=
\begin{bmatrix}
-1 & -1 & \cdots &1 \\
1 & 0 & \cdots & 1 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{bmatrix}
```

- $V$ possui os mesmos autovetores $v_i$ de P

$A$ assim como $P$ não são invertíveis. $A$ é uma matriz de projeção, então intuitivamente não é lógico dizer que A possui inversa. $P$ de mesmo modo, projeta vetores $v$ no espaço ortogonal a **1**, além disso, $P$ possui autovalores nulos, o que indica que existem vetores $v \neq 0$ tais que $Pv = 0$ (portanto $P$ não é injetiva, e consequentemente não invertível).

Com uma analogia dos hexágonos, se $P$ fosse invertível, seria possível para qualquer hexágono distante das proximidades da origem ser transladado para a origem e depois retornar ao mesmo local, porém isso é impossível porque existem infinitos locais que o hexágono poderia estar antes de ser aproximado à origem.

Para os vetores $v$ que não possuem componentes em **1**, a projeção não modificará nada no vetor. Da mesma maneira que, dados dois vetores $u$ e $v$ já no espaço ortogonal a **1**, não haverá distorções em suas distâncias relativas. Como $P$ diminui para cada vetor a média de seus valores, os vetores apenas serão transladados.

