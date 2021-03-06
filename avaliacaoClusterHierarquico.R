###############################################################################
# Author: Felipe Albuquerque
###############################################################################

#
# Cria��o da classe que representa um k-grupo. Guarda uma matriz de dados no atributo "grupo"
#
setClass("k.grupo", representation(grupo = "matrix"))

#
# Valida o valor de k, que deve estar entre 1 e N. 
# Retorna 1 quando k � inv�lido e 0 quando � v�lido
#
isKValido = function(k, n) {((k > 1) & (k < n))}

#
# Divide os dados em k grupos
# 
# dados: os dados que devem ser divididos
# k: o n�mero de grupos. Deve estar entre 1 e o n�mero de dados
#
dividirEmKGrupos = function(dados, k) {
    
    # Valida��o de k
    linhas = length(dados[,1])
    
    if (!isKValido(k, linhas)) {
        stop(paste("O valor de k deve estar entre 1 e", linhas))
    }
    
    # Inicializa��o do retorno. Cria��o das estruturas dos grupos
    cols = length(dados[1,])
    conjuntoRetornado = c()
    
    for (j in 1:k) {
        conjuntoRetornado = c(conjuntoRetornado, new("k.grupo", grupo = matrix(nrow = 0, ncol = cols)))
    }
    
    # Cria��o dos grupos    
    
    for (i in 1:linhas) {
        pos = i %% k # Sele��o do grupo baseada no resto da divis�o
        
        if (pos == 0) { # Como o �ndice come�a em "1", utilizar "k" ao inv�s de "0"
            pos = k
        }
        
        # Adiciona o dado ao grupo selecionado
        conjuntoRetornado[[pos]]@grupo = rbind(conjuntoRetornado[[pos]]@grupo, dados[i,])
    }
    
    conjuntoRetornado
}

#
# Retorna a matriz sem a coluna de classe
# 
# m: a matriz
# indiceColunaClasse: o �ndice da coluna de classe
#
extrairMatrizSemColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    cols = length(m[1,])
    
    if (indiceColunaClasse == cols) { # coluna de classe � a �ltima
        m[1:linhas, (1:cols - 1)]
    } else {
        x = m[1:linhas, 1:(indiceColunaClasse -1)]
        cbind(m[1:linhas, (indiceColunaClasse + 1):cols])
    }
    
}

#
# Retorna a coluna de classe da matriz
# 
# m: a matriz
# indiceColunaClasse: o �ndice da coluna de classe
#
extrairColunaDeClasse = function(m, indiceColunaClasse) {
    linhas = length(m[,1])
    as.matrix(m[1:linhas, indiceColunaClasse:indiceColunaClasse], nrow = linhas, ncol = 1, byrow = TRUE)
}

#
# Avalia o classificador baseado em agrupamento hier�rquico utilizando a t�cnica de cross-validation
#
avaliarAgrupamentoHierarquico = function(dados, k, indiceColunaClasse, funcaoMedicaoDistancia) {
    m = as.matrix(dados)
    linhas = length(m[,1])
    cols = length(m[1,])
    
    # Valida��o de k
    if (!isKValido(k, linhas)) {
        stop(paste("O valor de k deve estar entre 1 e", linhas))
    }
    
    # Execu��o da avalia��o
    grupos = dividirEmKGrupos(m, k)
    
    print(paste("Valida��o cruzada para k =", k))
    print(paste("Fun��o de medi��o de dist�ncia selecionada:", funcaoMedicaoDistancia))
    
    tp = 0
    fp = 0
    fn = 0
    tn = 0
    total = 0
    
    for (i in 1:k) {
        
        # Escolhe o grupo de testes
        grupoTeste = i
        conjuntoTestes = grupos[[grupoTeste]]@grupo
        
        # Como � n�o supervisionado, n�o h� modelo de treinamento
    
        # Teste
        x = extrairMatrizSemColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        y = extrairColunaDeClasse(conjuntoTestes, indiceColunaClasse)
        
        cluster = hclust(dist(x), funcaoMedicaoDistancia)
        
        # Algoritmo preparado apenas para classifica��o bin�ria
        resultadoTestes = table(cutree(cluster, k = 2), y) 
        
        # Incremento dos valores de medidas
        tp = tp + resultadoTestes[1, 1]
        fp = fp + resultadoTestes[1, 2]
        fn = fn + resultadoTestes[2, 1]
        tn = tn + resultadoTestes[2, 2]
        total = total + length(conjuntoTestes[,1])
    }
    
    # C�lculos das medidas
    p = tp + fn
    n = tn + fp
    taxaFalsaAceitacao = fp / total
    taxaFalsaRejeicao = fn / total
    sensibilidade = tp / p
    especificidade = tn / n
    erroMedio = (fp + fn) / total
    
    print(paste("Positivos:", p))
    print(paste("Negativos:", n))
    print(paste("Verdadeiros Positivos:", tp))
    print(paste("Falsos Positivos:", fp))
    print(paste("Falsos Negativos:", fn))
    print(paste("Verdadeiros Negativos:", tn))
    print(paste("Erro m�dio:", erroMedio))
    print(paste("Taxa de Falsa Aceita��o:", taxaFalsaAceitacao))
    print(paste("Taxa de Falsa Rejei��o:", taxaFalsaRejeicao))
    print(paste("Sensibilidade:", sensibilidade))
    print(paste("Especificidade:", especificidade))
}

# Carga dos dados
f = file("D:\\USP\\PPgSI - Mestrado\\2010-2_Sem\\SIN 5007-1 - Reconhecimento de Padr�es\\Atividades\\Atividade 6\\census_income_colunas-5_6-selecionadas.csv")
censusIncome = read.table(f, sep = ",", colClasses = c("integer", "integer", "integer"), header = TRUE)

require("graphics")

# Testes para diferentes fun��es de medi��o de dist�ncia 
funcoes = c("ward", "single", "complete", "average", "mcquitty", "centroid")

for (func in funcoes) {
    avaliarAgrupamentoHierarquico(dados = censusIncome, k = 10, indiceColunaClasse = 3, funcaoMedicaoDistancia = func)
}

plot(wardCluster)
plot(singleCluster)
plot(completeCluster)
plot(averageCluster)
plot(mcquittyCluster)
plot(centroidCluster)