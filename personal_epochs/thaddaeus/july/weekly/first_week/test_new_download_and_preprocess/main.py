from sunnyhills.pipeline_functions import better_download, better_preprocess

tic_id = 'TIC_17399504'

current_dir = 'personal_epochs/thaddaeus/july/weekly/first_week/test_new_download_and_preprocess/'


_, string_one, string_two, string_three, last_dates = better_download(tic_id, current_dir)

print(string_one, string_two, string_three)

better_preprocess(tic_id, current_dir+tic_id+'.csv', save_directory=current_dir, last_dates=last_dates)