nohup python3 train_model.py --model_name model_full_K562 --model_dir ./ --profile_main_dir ../../model_PRIMEloci/1_prep_data_for_training/K562_profiles > model_full_K562.log 2>&1 &
nohup python3 train_model.py --model_name model_full_K562_wt10M --model_dir ./ --profile_main_dir ../../model_PRIMEloci/1_prep_data_for_training/K562_wt10M_profiles > model_full_K562_wt10M.log 2>&1 &
nohup python3 train_model.py --model_name model_full_GM12878 --model_dir ./ --profile_main_dir ../../model_PRIMEloci/1_prep_data_for_training/GM12878_profiles > model_full_GM12878.log 2>&1 &
nohup python3 train_model.py --model_name model_full_GM12878_wt10M --model_dir ./ --profile_main_dir ../../model_PRIMEloci/1_prep_data_for_training/GM12878_wt10M_profiles > model_full_GM12878_wt10M.log 2>&1 &
