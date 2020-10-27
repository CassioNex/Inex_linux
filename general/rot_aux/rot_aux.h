/*
* rot_aux.h
*
* Created: 05/08/2014 14:21:51
*  Author: User
*/

#ifndef _ROT_AUX_H_
#define _ROT_AUX_H_

#include <stdio.h>
#include <stdarg.h>

#include <asf.h>
#include "common_def.h"

//! [definição de constantes]
// Buffer ASCII:
#define TAM_BUFFER_ASC                      8

// Opções de formatação de string numérica:
#define DESPREZA_ZEROS_ESQUERDA             0
#define MANTEM_ZEROS_ESQUERDA               1

// Flags para as funções numero_str e numero_buffer:
#define FLAGS_NUMERO_CONSIDERA_SINAL        0x01
#define FLAGS_NUMERO_CONSIDERA_DP           0x02
//! [definição de constantes]

//! [definição de estruturas do módulo]
//! [definição de estruturas do módulo]

//! [Variáveis globais públicas em memória de dados - usar diretiva extern]
//! [flags]
extern PUBLIC union unsigned_char rotaux_flags;
/*
* bit0:  F_AUX
* bit1:
* bit2:
* bit3:
* bit4:
* bit5:
* bit6:
* bit7:
*/
#define F_AUX	rotaux_flags.bit0
//! [flags]

extern PUBLIC unsigned char buffer_ascii[ TAM_BUFFER_ASC ];
//! [Variáveis globais públicas em memória de dados - usar diretiva extern]

//! [Variáveis globais públicas em memória de programa - usar diretiva extern]
//! [Variáveis globais públicas em memória de programa - usar diretiva extern]

//! [macros]
#define total_elementos_array( a )          ( sizeof( a ) / sizeof( a[ 0 ] ) )
//! [macros]

//! [prototipagem de funções públicas]
PUBLIC void copia_bytes( unsigned char *origem, unsigned char *destino, unsigned int tamanho );
PUBLIC void copia_bytes_flash( const unsigned char *origem, unsigned char *destino, unsigned int tamanho );
PUBLIC void copia_string( unsigned char *origem, unsigned char *destino );
PUBLIC void copia_string_flash( const unsigned char *origem, unsigned char *destino );

PUBLIC void concatena_string( unsigned char *origem, unsigned char *destino );
PUBLIC void concatena_string_flash( const unsigned char *origem, unsigned char *destino );

PUBLIC unsigned char lower_chr( unsigned char c );
PUBLIC unsigned char upper_chr( unsigned char c );

PUBLIC void lower_str( unsigned char *str );
PUBLIC void upper_str( unsigned char *str );

PUBLIC int pos_caracter_buffer( unsigned char caracter, unsigned char *buffer, unsigned int total_bytes );
PUBLIC int pos_str_buffer( unsigned char *str, unsigned char *buffer, unsigned int total_bytes );
PUBLIC int pos_str_flash_buffer( const unsigned char *str, unsigned char *buffer, unsigned int total_bytes );
PUBLIC unsigned char compara_str( unsigned char *str1, unsigned char *str2 );
PUBLIC unsigned char compara_str_flash( const unsigned char *str1, unsigned char *str2 );
PUBLIC unsigned char numero_str( unsigned char *str, unsigned char flags );
PUBLIC unsigned char numero_buffer( unsigned char *buffer, unsigned char total, unsigned char flags );
PUBLIC void preenche_buffer( unsigned char *buffer, unsigned char valor, unsigned int qtde );
PUBLIC unsigned char compara_buffer( unsigned char *buffer_a, unsigned char *buffer_b, unsigned int qtde );
PUBLIC unsigned char compara_flash_buffer( const unsigned char *buffer_a, unsigned char *buffer_b, unsigned int qtde );

PUBLIC unsigned int tamanho_str( unsigned char *str );
PUBLIC unsigned int tamanho_str_flash( const unsigned char *str );
PUBLIC void formata_str_valor( long numero, unsigned char casas_decimais, unsigned char max_caracteres, unsigned char opcao, unsigned char *saida );
PUBLIC unsigned int formata_str( unsigned char *saida, const unsigned char *str, ... );

PUBLIC inline unsigned char retorna_byte_high( unsigned int w );
PUBLIC inline unsigned char retorna_byte_low( unsigned int w );
PUBLIC inline unsigned int retorna_word( unsigned char h, unsigned char l );

PUBLIC long arredonda_valor( float v );

PUBLIC unsigned long ascii_para_dec( unsigned char *ascii, unsigned char tamanho );
PUBLIC void dec_para_ascii( unsigned long dec, unsigned char *ascii );
PUBLIC unsigned char bcd_para_dec( unsigned char bcd );
PUBLIC unsigned char dec_para_bcd( unsigned char dec );

PUBLIC void dword_para_hex( unsigned long valor, unsigned char *hex );
PUBLIC unsigned long hex_para_dword( unsigned char *hex, unsigned char qtde_digitos );
PUBLIC void byte_para_hex( unsigned char valor, unsigned char *hex );
PUBLIC unsigned char hex_para_byte( unsigned char *hex );
//! [prototipagem de funções públicas]

#endif // _ROT_AUX_H_
