/*
* common_def.h
*/

#ifndef COMMON_DEF_H_
#define COMMON_DEF_H_

#include "var.h"

#define NAO     		0
#define	SIM         	1

#define PUBLIC
#define PRIVATE         static

//! [Tipo de app]
#define SEM_BIOMETRIA	0
#define	COM_BIOMETRIA  	1
#define TIPO_APP        COM_BIOMETRIA
//! [Tipo de app]

//! [Versao do algoritmo de biometria]
// V1_0: versão arm original
// V2_0: versão arm "step2-finish"
// V_PC: versão pc original
#define ALG_V1_0        0
#define ALG_V2_0        1
#define ALG_V_PC        2
#define VERSAO_ALG      ALG_V_PC
//! [Versao do algoritmo de biometria]

#endif /* COMMON_DEF_H_ */
