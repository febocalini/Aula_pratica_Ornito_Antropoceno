##Stairway Xipho brejos
### fazer o download do programa
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git

### run_stairway plot --
cd /home/fernanda/Downloads/stairway_plot_v2.1.2
#pop PCE
java -cp stairway_plot_es Stairbuilder /media/fernanda/KINGSTON/stairway_linux/PCE_xipho/PCE_xipho.blueprint
##This command will create within seconds new directories and a bash script.

##Next execute the bash script:
bash /media/fernanda/KINGSTON/stairway_linux/PCE_xipho/PCE_xipho.blueprint.sh # note: don’t forget the extension (.sh)

#pop CE
java -cp stairway_plot_es Stairbuilder /media/fernanda/KINGSTON/stairway_linux/Ceara_xipho/CE_xipho.blueprint
##This command will create within seconds new directories and a bash script.
##Next execute the bash script:
bash /media/fernanda/KINGSTON/stairway_linux/Ceara_xipho/CE_xipho.blueprint.sh # note: don’t forget the extension (.sh)


#pop AM
java -cp stairway_plot_es Stairbuilder /media/fernanda/KINGSTON/stairway_linux/AM_xipho/AM_xipho.blueprint

##This command will create within seconds new directories and a bash script.
##Next execute the bash script:
bash /media/fernanda/KINGSTON/stairway_linux/AM_xipho/AM_xipho.blueprint.sh # note: don’t forget the extension (.sh)

#pop AFS
java -cp stairway_plot_es Stairbuilder /media/fernanda/KINGSTON/stairway_linux/AF_xipho/AF_xipho.blueprint

##This command will create within seconds new directories and a bash script.
##Next execute the bash script:
bash /media/fernanda/KINGSTON/stairway_linux/AF_xipho.blueprint.sh # note: don’t forget the extension (.sh)
