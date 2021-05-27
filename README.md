# inovanex - Linux / h7 libs
Projetos de fw/sw - Instituto Inovanex - Linux dk2 board / libs h7 mcu

Pastas:

  Linux (testado na mpu stm32mp157c-dk2):
  
    - App_bio_demo_x --> aplicacoes em C++ (QT), para testes de biometria com sensor fpc1011f3, e com o algoritmo Chines;
    - Gpio_demo --> aplicacao em C para teste de botoes e leds no linux (interessante desse sw eh o uso de multithread posix);
    - Spi_demo  --> aplicacao em C++/C para teste de spi master em linux;
    - Libs --> pasta de bibliotecas
      - fpc1011f3 --> driver para o sensor fpc1011f3;
      - finger    --> algoritmo de biometria, em tres versoes: v_pc, v1_0 (arm), e v2_0 (arm). As funcionais sao a v_pc e a v1_0;
      - common    --> arquivos comuns a todas as libs;
      - resources --> aquivos de recursos, como imagens etc.

  H7 (família de mcus de alta performance):
  
    - Libs --> pasta de bibliotecas
      - fpc1011f3 --> driver para o sensor fpc1011f3;
      - finger    --> algoritmo de biometria, em tres versoes: v_pc, v1_0 (arm), e v2_0 (arm). As funcionais sao a v_pc e a v1_0;
      - common    --> arquivos comuns a todas as libs;
      - resources --> aquivos de recursos, como imagens etc.
