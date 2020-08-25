#ifndef VAR_H_
#define VAR_H_

//****************************************
union unsigned_char
{
	struct //char_1
	{
		unsigned char bit0:1;
		unsigned char bit1:1;
		unsigned char bit2:1;
		unsigned char bit3:1;
		unsigned char bit4:1;
		unsigned char bit5:1;
		unsigned char bit6:1;
		unsigned char bit7:1;
	};
	//-----
	struct //char_2
	{
		unsigned char bit1_0:2;
		unsigned char bit3_2:2;
		unsigned char bit5_4:2;
		unsigned char bit7_6:2;
	};
	//-----
	struct //char_3
	{
		unsigned char bit3_0:4;
		unsigned char bit7_4:4;
	};
	//-----
	struct //char_4
	{
		unsigned char nibble0:4;
		unsigned char nibble1:4;
	};
	//-----
	struct //char_5
	{
		unsigned char nibble_low:4;
		unsigned char nibble_high:4;
	};
	//-----
	struct //char_6
	{
		unsigned char value;
	};
	//-----
	struct //char_7
	{
		unsigned char uc;
	};
	//-----
	struct //char_8
	{
		unsigned char byte;
	};
	struct
	{
		unsigned char bits;
	};
};

//*********************************
union unsigned_int
{
	struct //int_1
	{
		unsigned char bit0:1;
		unsigned char bit1:1;
		unsigned char bit2:1;
		unsigned char bit3:1;
		unsigned char bit4:1;
		unsigned char bit5:1;
		unsigned char bit6:1;
		unsigned char bit7:1;
		unsigned char bit8:1;
		unsigned char bit9:1;
		unsigned char bit10:1;
		unsigned char bit11:1;
		unsigned char bit12:1;
		unsigned char bit13:1;
		unsigned char bit14:1;
		unsigned char bit15:1;
	};
	//-----
	struct //int_2
	{
		unsigned char bit1_0:2;
		unsigned char bit3_2:2;
		unsigned char bit5_4:2;
		unsigned char bit7_6:2;
		unsigned char bit9_8:2;
		unsigned char bit11_10:2;
		unsigned char bit13_12:2;
		unsigned char bit15_14:2;
	};
	//-----
	struct //int_3
	{
		unsigned char bit3_0:4;
		unsigned char bit7_4:4;
		unsigned char bit11_8:4;
		unsigned char bit15_12:4;
	};
	//-----
	struct //int_4
	{
		unsigned char nibble0:4;
		unsigned char nibble1:4;
		unsigned char nibble2:4;
		unsigned char nibble3:4;
	};
	//-----
	struct //int_9
	{
		unsigned char uc0;
		unsigned char uc1;
	};
	//-----
	struct //int_10
	{
		unsigned char uc_low;
		unsigned char uc_high;
	};
	//-----
	struct //int_11
	{
		unsigned char uc[2];
	};
	//-----
	struct //int_12
	{
		unsigned char byte0;
		unsigned char byte1;
	};
	//-----
	struct //int_13
	{
		unsigned char byte_low;
		unsigned char byte_high;
	};
	//-----
	struct //int_14
	{
		unsigned char byte[2];
	};
	//-----
	struct //int_15
	{
		unsigned int value;
	};
	//-----
	struct //int_16
	{
		unsigned int ui;
	};
	//-----
	struct //int_17
	{
		unsigned int word;
	};
};

//*********************************
union unsigned_long
{
	struct //long_1
	{
		unsigned char bit0:1;
		unsigned char bit1:1;
		unsigned char bit2:1;
		unsigned char bit3:1;
		unsigned char bit4:1;
		unsigned char bit5:1;
		unsigned char bit6:1;
		unsigned char bit7:1;
		unsigned char bit8:1;
		unsigned char bit9:1;
		unsigned char bit10:1;
		unsigned char bit11:1;
		unsigned char bit12:1;
		unsigned char bit13:1;
		unsigned char bit14:1;
		unsigned char bit15:1;
		unsigned char bit16:1;
		unsigned char bit17:1;
		unsigned char bit18:1;
		unsigned char bit19:1;
		unsigned char bit20:1;
		unsigned char bit21:1;
		unsigned char bit22:1;
		unsigned char bit23:1;
		unsigned char bit24:1;
		unsigned char bit25:1;
		unsigned char bit26:1;
		unsigned char bit27:1;
		unsigned char bit28:1;
		unsigned char bit29:1;
		unsigned char bit30:1;
		unsigned char bit31:1;
	};
	//-----
	struct //long_2
	{
		unsigned char bit1_0:2;
		unsigned char bit3_2:2;
		unsigned char bit5_4:2;
		unsigned char bit7_6:2;
		unsigned char bit9_8:2;
		unsigned char bit11_10:2;
		unsigned char bit13_12:2;
		unsigned char bit15_14:2;
		unsigned char bit17_16:2;
		unsigned char bit19_18:2;
		unsigned char bit21_20:2;
		unsigned char bit23_22:2;
		unsigned char bit25_24:2;
		unsigned char bit27_26:2;
		unsigned char bit29_28:2;
		unsigned char bit31_30:2;
	};
	//-----
	struct //long_3
	{
		unsigned char bit3_0:4;
		unsigned char bit7_4:4;
		unsigned char bit11_8:4;
		unsigned char bit15_12:4;
		unsigned char bit19_16:4;
		unsigned char bit23_20:4;
		unsigned char bit27_24:4;
		unsigned char bit31_28:4;
	};
	//-----
	struct //long_4
	{
		unsigned char nibble0:4;
		unsigned char nibble1:4;
		unsigned char nibble2:4;
		unsigned char nibble3:4;
		unsigned char nibble4:4;
		unsigned char nibble5:4;
		unsigned char nibble6:4;
		unsigned char nibble7:4;
	};
	//-----
	struct //long_9
	{
		unsigned char uc0;
		unsigned char uc1;
		unsigned char uc2;
		unsigned char uc3;
	};
	//-----
	struct //long_10
	{
		unsigned char uc_low;
		unsigned char uc_high;
		unsigned char uc_medium;
		unsigned char uc_upper;
	};
	//-----
	struct //long_11
	{
		
	};
	//-----
	struct //long_12
	{
		unsigned char byte0;
		unsigned char byte1;
		unsigned char byte2;
		unsigned char byte3;
	};
	//-----
	struct //long_13
	{
		unsigned char byte_low;
		unsigned char byte_high;
		unsigned char byte_medium;
		unsigned char byte_upper;
	};
	//-----
	struct //long_14
	{
		unsigned char byte[4];
	};
	//-----
	struct //long_18
	{
		unsigned int ui0;
		unsigned int ui1;
	};
	//-----
	struct //long_19
	{
		unsigned int ui_low;
		unsigned int ui_high;
	};
	//-----
	struct //long_20
	{
		unsigned int ui[2];
	};
	//-----
	struct //long_21
	{
		unsigned int word0;
		unsigned int word1;
	};
	//-----
	struct //long_22
	{
		unsigned int word_low;
		unsigned int word_high;
	};
	//-----
	struct //long_23
	{
		unsigned int word[2];
	};
	//-----
	struct //long_24
	{
		unsigned long value;
	};
	//-----
	struct //long_25
	{
		unsigned long ul;
	};
	//-----
	struct //long_26
	{
		unsigned long dword;
	};
};

#endif /* VAR_H_ */
