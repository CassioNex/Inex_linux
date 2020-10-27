/*
* rot_aux.c
*
* Created: 05/08/2014 14:21:51
*  Author: User
*/

#include "rot_aux.h"

//! [prototipagem de fun��es privadas ao m�dulo]
//! [prototipagem de fun��es privadas ao m�dulo]


//! [Vari�veis globais p�blicas em mem�ria de dados]
PUBLIC unsigned char buffer_ascii[ TAM_BUFFER_ASC ];
PUBLIC union unsigned_char rotaux_flags;
//! [Vari�veis globais p�blicas em mem�ria de dados]

//! [Vari�veis globais p�blicas em mem�ria de programa]
//! [Vari�veis globais p�blicas em mem�ria de programa]

//! [Vari�veis globais privadas em mem�ria de dados]
//! [Vari�veis globais privadas em mem�ria de dados]

//! [Vari�veis globais privadas em mem�ria de programa]
//! [Vari�veis globais privadas em mem�ria de programa]

//! [Fun��es p�blicas]
/**
 * \external
 * copia_bytes
 *
 * Rotina para copiar byte a byte os dados de um buffer para outro
 *
 * \param[in]	string de origem, quantidade de bytes
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void copia_bytes( unsigned char *origem, unsigned char *destino, unsigned int tamanho )
{
    while( tamanho-- )
        *destino++ = *origem++;
}

/**
 * \external
 * copia_bytes_flash
 *
 * Rotina para copiar byte a byte os dados de um buffer em ROM para outro em RAM
 *
 * \param[in]	string de origem, quantidade de bytes
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void copia_bytes_flash( const unsigned char *origem, unsigned char *destino, unsigned int tamanho )
{
    while( tamanho-- )
        *destino++ = *origem++;
}

/**
 * \external
 * copia_string
 *
 * Rotina para copiar byte a byte de uma string em RAM para um buffer em RAM
 *
 * \param[in]	string de origem
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void copia_string( unsigned char *origem, unsigned char *destino )
{
    while( *origem )
        *destino++ = *origem++;
    *destino = '\0';
}

/**
 * \external
 * copia_string_flash
 *
 * Rotina para copiar byte a byte de uma string em ROM para um buffer em RAM
 *
 * \param[in]	string de origem
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void copia_string_flash( const unsigned char *origem, unsigned char *destino )
{
    while( *origem )
        *destino++ = *origem++;
    *destino = '\0';
}

/**
 * \external
 * concatena_string
 *
 * Rotina para concatenar duas strings. O resultado � armazenado na string de destino
 *
 * \param[in]	string de origem
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void concatena_string( unsigned char *origem, unsigned char *destino )
{
    unsigned int i;

    i = tamanho_str( destino );
    while( *origem )
        destino[ i++ ] = *origem++;
    destino[ i ] = '\0';
}

/**
 * \external
 * concatena_string_flash
 *
 * Rotina para concatenar duas strings. O resultado � armazenado na string de destino
 *
 * \param[in]	string de origem
 * \param[out]	string de destino
 *
 * \return		void
 */
PUBLIC void concatena_string_flash( const unsigned char *origem, unsigned char *destino )
{
    unsigned int i;

    i = tamanho_str( destino );
    while( *origem )
        destino[ i++ ] = *origem++;
    destino[ i ] = '\0';
}

/**
 * \external
 * lower_chr
 *
 * Rotina para converter um caractere para min�sculo
 *
 * \param[in]	caractere
 * \param[out]	void
 *
 * \return		caractere em min�sculo (se estiver entre 'A' e 'Z'), sen�o o pr�prio caractere
 */
PUBLIC unsigned char lower_chr( unsigned char c )
{
    if( ( c >= 'A' ) && ( c <= 'Z' ) )
        return( c - 'A' + 'a'  );

    // 0xE0: '�', 0xC0: '�'
    // 0xDC: '�'
    if( ( c >= 0xC0 ) && ( c <= 0xDC ) )
        return( c - 0xC0 + 0xE0 );

    return( c );
}

/**
 * \external
 * upper_chr
 *
 * Rotina para converter um caractere para mai�sculo
 *
 * \param[in]	caractere
 * \param[out]	void
 *
 * \return		caractere em mai�sculo (se estiver entre 'a' e 'z'), sen�o o pr�prio caractere
 */
PUBLIC unsigned char upper_chr( unsigned char c )
{
    if( ( c >= 'a' ) && ( c <= 'z' ) )
        return( c - 'a' + 'A' );

    // 0xE0: '�', 0xC0: '�'
    // 0xFC: '�'
    if( ( c >= 0xE0 ) && ( c <= 0xFC ) )
        return( c - 0xE0 + 0xC0 );

    return( c );
}

/**
 * \external
 * lower_str
 *
 * Rotina para converter uma string para min�sculo
 *
 * \param[in]	string
 * \param[out]	string
 *
 * \return		void
 */
PUBLIC void lower_str( unsigned char *str )
{
    unsigned char c;

    while( ( c = *str ) )
    {
        *str = lower_chr( c );
        str++;
    }
}

/**
 * \external
 * upper_str
 *
 * Rotina para converter uma string para mai�sculo
 *
 * \param[in]	string
 * \param[out]	string
 *
 * \return		void
 */
PUBLIC void upper_str( unsigned char *str )
{
    unsigned char c;

    while( ( c = *str ) )
    {
        *str = upper_chr( c );
        str++;
    }
}

/**
 * \external
 * pos_caracter_buffer
 *
 * Rotina para retornar a posi��o inicial de um buffer onde se encontra determinado ca-
 * ractere
 *
 * \param[in]	caractere a ser localizado, buffer com poss�vel caractere e tamanho total
				do buffer
 * \param[out]	void
 *
 * \return		retorna a posi��o no buffer referente a posi��o do caractere localizado
 *				ou -1 caso n�o tenha encontrado o caractere
 */
PUBLIC int pos_caracter_buffer( unsigned char caracter, unsigned char *buffer, unsigned int total_bytes )
{
    unsigned int i;
    
    for( i = 0 ; i < total_bytes ; i++ )
    {
        if( buffer[ i ] == caracter )
            return( i );
    }
    
    return( -1 );
}

/**
 * \external
 * pos_str_buffer
 *
 * Rotina para retornar a posi��o inicial de um buffer onde se encontra determinada string
 *
 * \param[in]	string a ser localizada, buffer com poss�vel string e tamanho total do buffer
 * \param[out]	void
 *
 * \return		retorna a posi��o no buffer referente ao in�cio da string localizada ou
 *				-1 caso n�o tenha encontrado a string 
 */
PUBLIC int pos_str_buffer( unsigned char *str, unsigned char *buffer, unsigned int total_bytes )
{
    unsigned char aux;
    int           pos;
    unsigned int  i;
    unsigned int  j;
    unsigned char *ptr;
    
    aux = 0;
    i   = 0;
    while( i < total_bytes )
    {
        j   = 0;
        pos = -1;
        ptr = str;
        while( *ptr )
        {
            aux = 1;
            if( *ptr != buffer[ i + j ] )
            {
                aux = 0;
                break;
            }
            else
            {
                if( pos == -1 )
                {
                    pos = i;
                }
            }

            j++;
            ptr++;
            if( ( i + j ) > total_bytes )
            {
                return( -1 );
            }
        }
        
        if( aux )
        {
            return( pos );
        }
        
        i++;
    }
    
    return( -1 );
}

/**
 * \external
 * pos_str_flash_buffer
 *
 * Rotina para retornar a posi��o inicial de um buffer onde se encontra determinada string
 *
 * \param[in]	string a ser localizada, buffer com poss�vel string e tamanho total do buffer
 * \param[out]	void
 *
 * \return		retorna a posi��o no buffer referente ao in�cio da string localizada ou
 *				-1 caso n�o tenha encontrado a string 
 */
PUBLIC int pos_str_flash_buffer( const unsigned char *str, unsigned char *buffer, unsigned int total_bytes )
{
    unsigned char aux;
    int           pos;
    unsigned int  i;
    unsigned int  j;
    const unsigned char *ptr;
    
    aux = 0;
    i   = 0;
    while( i < total_bytes )
    {
        j   = 0;
        pos = -1;
        ptr = str;
        while( *ptr )
        {
            aux = 1;
            if( *ptr != buffer[ i + j ] )
            {
                aux = 0;
                break;
            }
            else
            {
                if( pos == -1 )
                {
                    pos = i;
                }
            }

            j++;
            ptr++;
            if( ( i + j ) > total_bytes )
            {
                return( -1 );
            }
        }
        
        if( aux )
        {
            return( pos );
        }
        
        i++;
    }
    
    return( -1 );
}

/**
 * \external
 * compara_str
 *
 * Rotina para comparar duas strings
 *
 * \param[in]	string 1 e string 2
 * \param[out]	void
 *
 * \return		retorna 1 se as strings forem id�nticas, 0 caso contr�rio
 */
PUBLIC unsigned char compara_str( unsigned char *str1, unsigned char *str2 )
{
    unsigned int tam;
    
    tam = tamanho_str( str1 );
    if( tamanho_str( str2 ) == tam )
    {
        if( pos_str_buffer( str1, str2, tam ) == 0 )
        {
            return( 1 );
        }
    }
    
    return( 0 );
}

/**
 * \external
 * compara_str_flash
 *
 * Rotina para comparar duas strings em mem�ria de programa
 *
 * \param[in]	string 1 e string 2
 * \param[out]	void
 *
 * \return		retorna 1 se as strings forem id�nticas, 0 caso contr�rio
 */
PUBLIC unsigned char compara_str_flash( const unsigned char *str1, unsigned char *str2 )
{
    unsigned int tam;
    
    tam = tamanho_str_flash( str1 );
    if( tamanho_str_flash( str2 ) == tam )
    {
        if( pos_str_flash_buffer( str1, str2, tam ) == 0 )
        {
            return( 1 );
        }
    }
    
    return( 0 );
}

/**
 * \external
 * numero_str
 *
 * Rotina para informar se a string passada � um n�mero ou n�o
 *
 * \param[in]	string a ser testada e flags de teste (FLAGS_NUMERO_STR_CONSIDERA_SINAL:
 *				deve ou n�o considerar sinal negativo na string / FLAGS_NUMERO_STR_CONSI-
 *				DERA_DP: deve ou n�o considerar ponto decimal na string)
 * \param[out]	void
 *
 * \return		retorna 1 se a string for um n�mero, 0 caso contr�rio
 */
PUBLIC unsigned char numero_str( unsigned char *str, unsigned char flags )
{
    unsigned char c;
    unsigned char i;
    unsigned char dp;
    
    i = 0;
    dp = 0;
    while( ( c = str[ i++ ] ) != '\0' )
    {
        if( ( c < '0' ) || ( c > '9' ) )
        {
            if( c == '-' )
            {
                if( ( flags & FLAGS_NUMERO_CONSIDERA_SINAL ) && ( i == 1 ) && ( str[ i ] != '\0' ) )
                {
                    continue;
                }
            }
            else if( ( c == '.' ) || ( c == ',' ) )
            {
                if( ( flags & FLAGS_NUMERO_CONSIDERA_DP ) && ( !dp ) )
                {
                    dp = 1;
                    continue;
                }
            }

            return( 0 );
        }
    }
    
    return( 1 );
}

/**
 * \external
 * numero_buffer
 *
 * Rotina para informar se o buffer passado � um n�mero ou n�o
 *
 * \param[in]	buffer a ser testado e total de bytes presentes
 * \param[out]	void
 *
 * \return		retorna 1 se o buffer for um n�mero, 0 caso contr�rio
 */
PUBLIC unsigned char numero_buffer( unsigned char *buffer, unsigned char total, unsigned char flags )
{
    unsigned char c;
    unsigned char i;
    unsigned char dp;
    
    i = 0;
    dp = 0;
    while( ( total-- ) && ( c = buffer[ i++ ] ) != '\0' )
    {
        if( ( c < '0' ) || ( c > '9' ) )
        {
            if( c == '-' )
            {
                if( ( flags & FLAGS_NUMERO_CONSIDERA_SINAL ) && ( i == 1 ) && ( total > 0 ) )
                {
                    continue;
                }
            }
            else if( ( c == '.' ) || ( c == ',' ) )
            {
                if( ( flags & FLAGS_NUMERO_CONSIDERA_DP ) && ( !dp ) )
                {
                    dp = 1;
                    continue;
                }
            }

            return( 0 );
        }
    }
    
    return( 1 );
}

/**
 * \external
 * preenche_buffer
 *
 * Rotina para preencher um buffer com determinado valor numa certa quantidade
 *
 * \param[in]	valor de preenchimento e quantidade de bytes que ser�o preenchidos
 * \param[out]	buffer que ser� preenchido
 *
 * \return		void
 */
PUBLIC void preenche_buffer( unsigned char *buffer, unsigned char valor, unsigned int qtde )
{
    while( qtde-- )
        *buffer++ = valor;
}

/**
 * \external
 * compara_buffer
 *
 * Rotina para comparar dois buffers e retornar o resultado 0 para diferente e 1 para
 * igual
 *
 * \param[in]	buffers a serem comparados e tamanho
 * \param[out]	void
 *
 * \return		0 para buffers diferentes e 1 para buffers iguais
 */
PUBLIC unsigned char compara_buffer( unsigned char *buffer_a, unsigned char *buffer_b, unsigned int qtde )
{
    while( qtde-- )
    {
        if( *buffer_a++ != *buffer_b++  )
        {
            return( 0 );
        }
    }

    return( 1 );
}

/**
 * \external
 * compara_flash_buffer
 *
 * Rotina para comparar dois buffers e retornar o resultado 0 para diferente e 1 para
 * igual
 *
 * \param[in]	buffers a serem comparados e tamanho
 * \param[out]	void
 *
 * \return		0 para buffers diferentes e 1 para buffers iguais
 */
PUBLIC unsigned char compara_flash_buffer( const unsigned char *buffer_a, unsigned char *buffer_b, unsigned int qtde )
{
    while( qtde-- )
    {
        if( *buffer_a++ != *buffer_b++  )
        {
            return( 0 );
        }
    }

    return( 1 );
}

/**
 * \external
 * tamanho_str
 *
 * Rotina para calcular o tamanho de uma string em mem�ria de dados
 *
 * \param[in]	string
 * \param[out]	void
 *
 * \return		n�mero de caracteres da string
 */
PUBLIC unsigned int tamanho_str( unsigned char *str )
{
    unsigned char tamanho = 0;
    
    while( str[ tamanho ] )
        tamanho++;
        
    return( tamanho );
}

/**
 * \external
 * tamanho_str_flash
 *
 * Rotina para calcular o tamanho de uma string em mem�ria de programa
 *
 * \param[in]	string
 * \param[out]	void
 *
 * \return		n�mero de caracteres da string
 */
PUBLIC unsigned int tamanho_str_flash( const unsigned char *str )
{
    unsigned char tamanho = 0;
    
    while( str[ tamanho ] )
        tamanho++;
        
    return( tamanho );
}

/**
 * \external
 * formata_str_valor
 *
 * Rotina para formatar um n�mero em uma string
 *
 * \param[in]	n�mero a escrever, n�mero de casas decimais, n�mero de colunas que cons-
 *				tituem a �rea de escrita (conta os d�gitos a partir dos menos significa-
 *				tivos, isto �, da direita para a esquerda) ou se igual a zero ent�o deve
 *				exibir todos os d�gitos, op��o de desprezar ou manter zeros � esquerda
 * \param[out]	string de sa�da
 *
 * \return		void
 */
PUBLIC void formata_str_valor( long numero, unsigned char casas_decimais, unsigned char max_caracteres, unsigned char opcao, unsigned char *saida )
{
    unsigned char i;
    unsigned char j;
    unsigned char pos_fim;
    unsigned char digitos  = 0;
    unsigned char negativo = 0;
    
    static unsigned char buffer[ TAM_BUFFER_ASC + 3 ];// 1 caracter para sinal, 1 caracter para ponto decimal e outro para o fim da string
    
    // Se o n�mero for negativo, seta flag e trabalha com seu m�dulo
    if( numero < 0 )
    {
        numero  *= -1;
        negativo =  1;
    }
    
    dec_para_ascii( numero, &buffer[ 1 ] ); 
    
    // Se desejar casas decimais, pula dois bytes ap�s os d�gitos ASCII: um byte � para o sinal
    // negativo (caso necess�rio) e o outro � para o ponto decimal
    if( casas_decimais > 0 )
    {
        pos_fim = TAM_BUFFER_ASC + 2;               
    }
    else
    {
        pos_fim = TAM_BUFFER_ASC + 1;
    }
    buffer[ pos_fim ] = '\0';
    
    // Checagem de limite do par�metro max_caracteres
    if( ( max_caracteres == 0 ) || ( max_caracteres > pos_fim ) )
    {
        max_caracteres = pos_fim;
    }
    
    // Checagem de limite do par�metro casas_decimais
    if( casas_decimais > ( TAM_BUFFER_ASC - 1 ) )
    {
        casas_decimais = ( TAM_BUFFER_ASC - 1 );
    }
    
    // Contabiliza d�gitos: se for zero, incrementa uma vez. Se for diferente de zero, incrementa
    // quantas vezes � divis�vel por 10
    if( numero == 0 )
    {
        digitos++;
    }
    else
    {
        while( numero )
        {
            digitos++;
            numero /= 10;
        }
    }
    
    // Se a quantidade de d�gitos puros for menor que a de casas decimais, adiciona a diferen�a
    if( digitos <= casas_decimais )
    {
        digitos += ( casas_decimais - digitos + 1 );
    }
    
    // Contabiliza d�gitos: se tiver casa decimal, incrementa
    if( casas_decimais > 0 )
    {
        digitos++;
    }
    
    // Contabiliza d�gitos: se for negativo, incrementa
    if( negativo )
    {
        digitos++;
    }
    
    // Se for para manter os zeros a esquerda, considera os d�gitos menos significativos,
    // isto �, da direita para a esquerda. Consequentemente, pelo fato do �ltimo la�o ler
    // o buffer de pos_fim - digitos, digitos precisa iniciar com max_caracteres
    if( opcao == MANTEM_ZEROS_ESQUERDA )
    {
        digitos = max_caracteres;
        
        // Se o n�mero n�o for negativo, o primeiro �ndice de buffer, destinado para o sinal,
        // estar� com o valor 0x00 e se n�o for subtra�do 1 de digitos, o �ltimo la�o n�o ser�
        // realizado
        if( digitos == pos_fim )
        {
            if( !negativo )
            {
                digitos--;
            }
        }
    }
    
    // Insere o sinal '-' se o n�mero for negativo, na posi��o adequada
    if( negativo )
    {
        buffer[ pos_fim - digitos ] = '-';
    }
    
    // Se desejar casas decimais, ent�o desloca de uma posi��o para a direita todos os d�gitos
    // ap�s a casa decimal, a fim de inserir a casa decimal na posi��o correta
    if( casas_decimais > 0 )
    {
        i = ( TAM_BUFFER_ASC + 1 );
        while( casas_decimais-- )
        {
            buffer[ i ] = buffer[ i - 1 ];
            i--;
        }
        
        buffer[ i ] = '.';
    }
    
    // Come�a a string a partir do sinal de menos (se o valor for negativo) ou do primeiro d�gito
    // do n�mero (se for positivo)
    i = 0;
    while( ( j = buffer[ ( pos_fim - digitos ) + i++ ] ) && ( max_caracteres > 0 ) )
    {
        *saida++ = j;
        max_caracteres--;
    }
    *saida = '\0';
}

/**
 * \external
 * formata_str
 *
 * Rotina que realiza a formata��o de uma string, podendo concatenar, transformar n�me-
 * ros, inserir casas decimais, zeros � esquerda etc, conforme sprintf()
 *
 * \param[in]	string contendo a formata��o (podendo incluir %d, %ld, %f, %s, %u etc)
 *				e os par�metros correspondentes aos caracteres de formata��o
 * \param[out]	buffer de sa�da da string
 *
 * \return		total de bytes presentes na string de sa�da
 */
PUBLIC unsigned int formata_str( unsigned char *saida, const unsigned char *str, ... )
{
    unsigned int n;
    va_list ap;

    va_start( ap, str );
    n = vsprintf( ( char * )saida, ( const char * )str, ap );
    va_end( ap );

    return( n );
}

/**
 * \external
 * retorna_byte_high
 *
 * Rotina para retornar o byte mais significativo de uma word
 *
 * \param[in]	word
 * \param[out]	void
 *
 * \return		byte mais significativo
 */
PUBLIC inline unsigned char retorna_byte_high( unsigned int w )
{
    return( w >> 8 );
}

/**
 * \external
 * retorna_byte_low
 *
 * Rotina para retornar o byte menos significativo de uma word
 *
 * \param[in]	word
 * \param[out]	void
 *
 * \return		byte menos significativo
 */
PUBLIC inline unsigned char retorna_byte_low( unsigned int w )
{
    return( w & 0xFF );
}

/**
 * \external
 * retorna_word
 *
 * Rotina para retornar uma word a partir de dois bytes
 *
 * \param[in]	byte mais e menos significativo
 * \param[out]	void
 *
 * \return		word
 */
PUBLIC inline unsigned int retorna_word( unsigned char h, unsigned char l )
{
    return( ( ( unsigned int )h << 8 ) + l );
}

/**
 * \external
 * arredonda_valor
 *
 * Rotina para arredondar um valor em ponto flutuante para um valor em inteiro, anali-
 * sando a primeira casa decimal
 *
 * \param[in]	valor em ponto flutuante
 * \param[out]	void
 *
 * \return		valor inteiro, arredondado
 */

PUBLIC long arredonda_valor( float v )
{
    if( v < 0.0 ) {
        return( v - 0.5 );
	}
	return( v + 0.5 );
}

/**
 * \external
 * ascii_para_dec
 *
 * Rotina para converter d�gitos em ASCII para seu correspondente valor inteiro
 *
 * \param[in]	buffer com valor em ascii e n�mero de d�gitos
 * \param[out]	void
 *
 * \return		valor inteiro
 */
PUBLIC unsigned long ascii_para_dec( unsigned char *ascii, unsigned char tamanho )
{
    unsigned long dec = 0;
    
    while( tamanho-- )
    {
        dec *= 10;
        dec += ( *ascii++ ) - '0';
    }
    
    return( dec );
}

/**
 * \external
 * dec_para_ascii
 *
 * Rotina para converter um valor inteiro de no m�ximo 32bits para seus correspondentes
 * d�gitos em ASCII
 *
 * \param[in]	valor inteiro
 * \param[out]	pointer para buffer de sa�da em ascii
 *
 * \return		void
 */
PUBLIC void dec_para_ascii( unsigned long dec, unsigned char *ascii )
{
    ascii[ 0 ] = '0';
    ascii[ 1 ] = '0';
    ascii[ 2 ] = '0';
    ascii[ 3 ] = '0';
    ascii[ 4 ] = '0';
    ascii[ 5 ] = '0';
    ascii[ 6 ] = '0';
    ascii[ 7 ] = '0';
    
    if( !dec )
        return;
    
    while( dec >= 100000000 )
    {
        dec -= 100000000;
    }
    while( dec >= 10000000 )
    {
        ascii[ 0 ]++;
        dec -= 10000000;
    }
    while( dec >= 1000000 )
    {
        ascii[ 1 ]++;
        dec -= 1000000;
    }
    while( dec >= 100000 )
    {
        ascii[ 2 ]++;
        dec -= 100000;
    }
    while( dec >= 10000 )
    {
        ascii[ 3 ]++;
        dec -= 10000;
    }
    while( dec >= 1000 )
    {
        ascii[ 4 ]++;
        dec -= 1000;
    }
    while( dec >= 100 )
    {
        ascii[ 5 ]++;
        dec -= 100;
    }
    while( dec >= 10 )
    {
        ascii[ 6 ]++;
        dec -= 10;
    }
    ascii[ 7 ] += dec;
}

/**
 * \external
 * bcd_para_dec
 *
 * Rotina para convers�o de n�mero BCD em n�mero decimal
 *
 * \param[in]	n�mero BCD a converter
 * \param[out]	void
 *
 * \return		n�mero decimal
 */
PUBLIC unsigned char bcd_para_dec( unsigned char bcd )
{
    return( ( ( bcd / 0x10 ) * 10 ) + ( bcd % 0x10 ) );
}

/**
 * \external
 * converte_dec_bcd
 *
 * Rotina para convers�o de n�mero decimal em n�mero BCD
 *
 * \param[in]		n�mero decimal a converter
 * \param[out]	void
 *
 * \return		n�mero BCD
 */
PUBLIC unsigned char dec_para_bcd( unsigned char dec )
{
    return( ( ( dec / 10 ) * 0x10 ) + ( dec % 10 ) );
}

/**
 * \external
 * dword_para_hex
 *
 * Rotina de convers�o de n�mero de 4 bytes para 8 d�gitos hexa
 *
 * \param[in]	n�mero a converter
 * \param[out]	buffer de sa�da
 *
 * \return		void
 */
PUBLIC void dword_para_hex( unsigned long valor, unsigned char *hex )
{
    union unsigned_long bytes;
    
    bytes.value = valor;
    
    byte_para_hex( bytes.byte[ 3 ], &hex[ 0 ] );
    byte_para_hex( bytes.byte[ 2 ], &hex[ 2 ] );
    byte_para_hex( bytes.byte[ 1 ], &hex[ 4 ] );
    byte_para_hex( bytes.byte[ 0 ], &hex[ 6 ] );
}

/**
 * \external
 * hex_para_dword
 *
 * Rotina de convers�o de no m�ximo 8 d�gitos hexa para n�mero de 4 bytes
 *
 * \param[in]	buffer com os d�gitos hexa e a quantidade de d�gitos
 * \param[out]	void
 *
 * \return		n�mero de 4 bytes
 */
PUBLIC unsigned long hex_para_dword( unsigned char *hex, unsigned char qtde_digitos )
{
    unsigned char nibble;
    unsigned long resultado = 0;
    
    while( qtde_digitos-- )
    {
        if( *hex >= 'A' )
        {
            nibble = 0x0A + *hex - 'A';
        }
        else
        {
            nibble = *hex - '0';
        }
        
        resultado *= 0x10;
        resultado += nibble;
        
        hex++;
    }
    
    return( resultado );
}

/**
 * \external
 * byte_para_hex
 *
 * Rotina de convers�o de n�mero de 1 byte para 2 d�gitos hexa
 *
 * \param[in]	n�mero a converter
 * \param[out]	buffer com 2 d�gitos hexa
 *
 * \return		void
 */
PUBLIC void byte_para_hex( unsigned char valor, unsigned char *hex )
{
    union unsigned_char byte;
    
    byte.value = valor;
    
    if( byte.nibble_high >= 0x0A )
    {
        hex[ 0 ] = 'A' + byte.nibble_high - 0x0A;
    }
    else
    {
        hex[ 0 ] = '0' + byte.nibble_high;
    }
    
    if( byte.nibble_low >= 0x0A )
    {
        hex[ 1 ] = 'A' + byte.nibble_low - 0x0A;
    }
    else
    {
        hex[ 1 ] = '0' + byte.nibble_low;
    }
}

/**
 * \external
 * hex_para_byte
 *
 * Rotina de convers�o de 2 d�gitos hexa para n�mero de 1 byte
 *
 * \param[in]	buffer contendo os dois d�gitos em hexa
 * \param[out]	void
 *
 * \return		n�mero de 1 byte
 */
PUBLIC unsigned char hex_para_byte( unsigned char *hex )
{
    unsigned char i;
    unsigned char nibble;
    unsigned char resultado = 0;
    
    for( i = 0 ; i < 2 ; i++ )
    {
        if( hex[i] >= 'A' )
        {
            nibble = 0x0A + hex[ i ] - 'A';
        }
        else
        {
            nibble = hex[ i ] - '0';
        }
        
        resultado *= 0x10;
        resultado += nibble;
    }
    
    return( resultado );
}
//! [Fun��es p�blicas]

//! [Fun��es privadas]
//! [Fun��es privadas]
