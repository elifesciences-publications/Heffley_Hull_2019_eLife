%WF CRP list of days

animals = {'img86', [], 'img88', 'img95', 'img96', 'img98', 'img99'};

days_1 = {'170321_img86', [], '170408_img88', '170605_img95', '170612_img96', '170628_img98', '170705_img99'};  %'170408_img87' first day of imaging but animal quit licking part way through
days_1_ROI = {[1,5], [], [1,2,5], [1:4], [1,2,4], [1:5], [3,4,5]};
days_1_stable_int = {[1,150], [], [1126,1200], [1417,1485], [808,890], [273,324], [230,292]};

days_post = {'170325_img86', [], '170415_img88', '170614_img95' '170622_img96', '170704_img98', '170710_img99'}; %'170326_img86',  are also 10%omission post learning
days_post_ROI = {[1,5], [], [1,2,5], [1:4], [1,2,4], [1:5], [3,4,5]};

days_UR = {'170327_img86', [], '170417_img88', [], '170623_img96', '170705_img98', '170711_img99'}; %img86_10%OR_10%%UR  170420_img87_10%/10%   170417_img88_10/10   '170420_img87'
days_UR_ROI = {[1,5], [], [1,2,5], [], [1,2,4], [1:5], [3,4,5]};


%% 
days_to_learn_criteria =     [3, 6, 5, 10, 7, 7, 6, 4, 8, ...
    7, 3, 4, 5, 6, 7, 3, 10, 11, 7, ...  %consider changing img98 to day 4
    8, 3, 4, 4, 3, 4, 4, 5];   %listed in numerical order by animal ID
days_to_PL_imaging_session = [4, 7, 6, 11, 8, 8, 8, 5, 10, ...
    8, 5, 5, 6, 7, 8, 4, 13, 12, 8, ...
    9, 5, 8, 5, 5, 5, 5, 6];
all_bx_animals = {'img86', 'img88', 'img89', 'img90', 'img91', 'img92', 'img93', 'img94', 'img95', ...
    'img96', 'img98', 'img99', 'img044', 'img050', 'img067', 'img070', 'img077', 'img081', 'img084', ...
    'img085', 'img087', 'img088', 'img089', 'img091', 'img092', 'img093', 'img094'};
