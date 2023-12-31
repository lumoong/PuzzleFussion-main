"""
Train a diffusion model on images.
"""

import argparse
import wandb
import torch as th
from jigsawplan import logger, dist_util
from jigsawplan.crosscut_dataset import load_crosscut_data, load_crosscut2_data
from jigsawplan.resample import create_named_schedule_sampler
from jigsawplan.script_util import (
    model_and_diffusion_defaults,
    create_model_and_diffusion,
    args_to_dict,
    add_dict_to_argparser,
    update_arg_parser,
)
from jigsawplan.train_util import TrainLoop


def main():

    wandb.init(
        # set the wandb project where this run will be logged
        project="postecha",
        # track hyperparameters and run metadata
        config={
        "learning_rate": 1e-5,
        "dataset": "plank2",
        "optimizer": "adamW",
        }
    ) 

    args = create_argparser().parse_args()
    update_arg_parser(args)

    dist_util.setup_dist()
    logger.configure(dir=f'ckpts/{args.exp_name}')

    logger.log("creating model and diffusion...")
    model, diffusion = create_model_and_diffusion(
        **args_to_dict(args, model_and_diffusion_defaults().keys())
    )
    model.to(dist_util.dev())
    schedule_sampler = create_named_schedule_sampler(args.schedule_sampler, diffusion)
    #rotation=True,
    '''args.num_channels = 256
    num_coords = 2 if not args.rotation else 4
    if args.dataset=='crosscut':
        args.input_channels = num_coords #+ (2*8)
        args.condition_channels = 66
        args.out_channels = num_coords * 1
    else:
        assert False, "DATASET NOT FOUND"'''
    '''res = dict(
            dataset='',
            use_checkpoint=False,
            input_channels=0,
            condition_channels=0,
            out_channels=0,
            num_channels=128
            )'''
    logger.log("creating data loader...")
    if args.dataset=='crosscut':
        data = load_crosscut_data(
            batch_size=args.batch_size,
            set_name=args.set_name,
            rotation=args.rotation,
            use_image_features=args.use_image_features,
        )
    if args.dataset=='crosscut2':
        data = load_crosscut2_data(
            batch_size=args.batch_size,
            set_name=args.set_name,
            rotation=args.rotation,
            use_image_features=args.use_image_features,
        )
    else:
        print('dataset not exist!')
        assert False

    logger.log("training...")
    TrainLoop(
        model=model,
        diffusion=diffusion,
        data=data,
        batch_size=args.batch_size,
        microbatch=args.microbatch,
        lr=args.lr,
        ema_rate=args.ema_rate,
        log_interval=args.log_interval,
        save_interval=args.save_interval,
        resume_checkpoint=args.resume_checkpoint,
        use_fp16=args.use_fp16,
        fp16_scale_growth=args.fp16_scale_growth,
        schedule_sampler=schedule_sampler,
        weight_decay=args.weight_decay,
        lr_anneal_steps=args.lr_anneal_steps,
    ).run_loop()
    dist_util.cleanup()


def create_argparser():
    defaults = dict(
        dataset = '',
        schedule_sampler= "uniform", #"loss-second-moment", "uniform",
        lr=1e-3,
        weight_decay=0.0,
        lr_anneal_steps=0,
        batch_size=1,
        microbatch=-1,  # -1 disables microbatches
        ema_rate="0.9999",  # comma-separated list of EMA values
        log_interval=500,
        save_interval=25000,
        resume_checkpoint="",
        use_fp16=False,
        fp16_scale_growth=1e-3,
    )
    parser = argparse.ArgumentParser()
    defaults.update(model_and_diffusion_defaults())
    add_dict_to_argparser(parser, defaults)
    return parser


if __name__ == "__main__":
    main()
