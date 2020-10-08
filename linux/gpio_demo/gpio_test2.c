#include <stdio.h>
#include <unistd.h>
#include <gpiod.h>
 
int main(int argc, char *argv[])
{
    struct gpiod_chip *output_chip;
    struct gpiod_line *output_line;
    int line_value;
 
    /* open /dev/gpiochip0 */
    output_chip = gpiod_chip_open_by_number(0);
    
    /* find BLUE_LED pin */
    output_line = gpiod_chip_find_line(output_chip, "BLUE_LED");
 
    /* config as output and set a description */
    gpiod_line_request_output(output_line, "embarcadosDemo",
                  GPIOD_LINE_ACTIVE_STATE_HIGH);
 
    while (1) {
        line_value = !line_value;
        gpiod_line_set_value(output_line, line_value);
        sleep(1);
    }
    
    return 0;
}
 
