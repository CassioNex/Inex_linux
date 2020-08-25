#include <stdio.h>

void main()
{
    char *buffer;
    char buf[10];
    int file;
    
    file=spi_init("/dev/spidev0.0"); //dev
    
    buf[0] = 0x41;
    buf[1] = 0xFF;
    spi_write(0xE6,0x0E,2,buf,file); //this will write value 0x41FF to the address 0xE60E
    
    buffer=(char *)spi_read(0xE6,0x0E,4,file); //reading the address 0xE60E
    
    close(file);
}
