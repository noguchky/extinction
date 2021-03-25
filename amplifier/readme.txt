python rbcpshell_v3.py 192.168.10.16 4660
でRBCP Shellに入る

RBCP shellのコマンド
rd 0xXX n
read コマンド アドレス0xXXから n番目のアドレスまでの値を読み込み
nは10進でも16進でも良い

wrb 0xXX 0x01
write コマンド 0xXXに0x01を書き込み

load XXXX.txt
load コマンド XXXX.txtにコマンドを記述しておくと、上から読み込んでくれる。

--------以下txtファイルの中身サンプル---------------------------
rd  0x10000000 48
wrb 0x10000010 0x01
wrb 0x10000011 0x41
rd  0x10000011 1
wrb 0x10000020 0x55
wrb 0x10000021 0x55
wrb 0x10000022 0x00
rd  0x10000020 8
wrb 0x10000012 0x01
rd  0x10000020 8
