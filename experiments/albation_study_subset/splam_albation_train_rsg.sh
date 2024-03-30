# for rsb in 1 2 3
# do 
#     rsg=5
#     echo python splam_albation_train.py --rsg $rsg --rsb $rsb
#     python splam_albation_train.py --rsg $rsg --rsb $rsb
# done

for rsg in 1 2 3 4
do 
    rsb=4
    echo python splam_albation_train.py --rsg $rsg --rsb $rsb
    python splam_albation_train.py --rsg $rsg --rsb $rsb
done

# rsg=5
# rsb=4

# echo python splam_albation_train.py --rsg $rsg --rsb $rsb
# python splam_albation_train.py --rsg $rsg --rsb $rsb
