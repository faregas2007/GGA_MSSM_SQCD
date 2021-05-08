sbotfile=('sbot.02' 'sbot.02.5' 'sbot.03' 'sbot.03.5' 'sbot.04' 'sbot.04.5' 'sbot.05' 'sbot.05.5' 'sbot.06' 'sbot.06.5' 'sbot.07' 'sbot.07.5' 'sbot.08' 'sbot.08.5' 'sbot.09' 'sbot.09.5' 'sbot.10' 'sbot.40')
stopfile=('stop.02' 'stop.02.5' 'stop.03' 'stop.03.5' 'stop.04' 'stop.04.5' 'stop.05' 'stop.05.5' 'stop.06' 'stop.06.5' 'stop.07' 'stop.07.5' 'stop.08' 'stop.08.5' 'stop.09' 'stop.09.5' 'stop.10' 'stop.40')
mArange=('1000' '1200' '1400' '1600' '1800' '2000' '2200' '2400' '2500' '2600' '2700' '2800')

for filename in ${sbotfile[@]}
do
 python tanbinput.py ${filename}
done

for filename in ${stopfile[@]}
do
 python tanbinput.py ${filename}
done

for mA in ${mArange[@]}
do
 python mAinput.py stop.40 ${mA}
 python mAinput.py stop.05 ${mA}
 python mAinput.py sbot.40 ${mA}
 python mAinput.py sbot.05 ${mA}
done
