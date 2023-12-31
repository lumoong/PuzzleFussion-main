MODEL_FLAGS="--dataset crosscut2 --batch_size 32 --set_name train"
TRAIN_FLAGS="--lr 1e-3 --save_interval 5000 --weight_decay 0.05 --log_interval 500"
#SAMPLE_FLAGS="--batch_size 1024 --num_samples 1024 --set_name test"

mpiexec -n 1 python image_train.py $MODEL_FLAGS $TRAIN_FLAGS --exp_name pred_s_wp_c_wimage 