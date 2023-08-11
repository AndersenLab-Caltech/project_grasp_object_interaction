%

clc
clear all
close all

nbrCues = 2;
ChanceLevel = 20*ones(nbrCues,1) + (rand(nbrCues,1)*5 - 1/2);
NotChance = 80*ones(nbrCues,1) + (rand(nbrCues,1)*15 - 5);

CI = 5*ones(nbrCues,1); %var(ChanceLevel)*ones(3,1)*10;
%ColorsCues = utile.get_color_rgb_codes({'redballoon', 'Auditory'});
ColorsCues = {[[0 51 51]/255],[[255 128 0]/256]};

figure();   

errorbar(1:length(NotChance), NotChance, CI,'-s', 'LineStyle', 'none','MarkerSize',15,...
        'MarkerEdgeColor',ColorsCues{2},'MarkerFaceColor',ColorsCues{2}, 'Color', ColorsCues{2}) 


hold on
errorbar(1:length(ChanceLevel), ChanceLevel, CI,'-s', 'LineStyle', 'none','MarkerSize',15,...
    'MarkerEdgeColor',ColorsCues{1},'MarkerFaceColor',ColorsCues{1}, 'Color', ColorsCues{1}) 
hold on

    
chance_level = 1/5*100;


l = line([0.5 3.5 ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
legend('Grasps', 'Pseudowords', 'Chance Level', 'fontsize', 12)

ylim([ 0 100]);
xticks([1:3])
xticklabels({'Image', 'Auditory', 'Written'})
%xtickangle(-45)


%%

ChanceLevel = [21 19]';
NotChance = [85 71]';

figure();


subplot(2,2,1)
bar([ChanceLevel, ChanceLevel + rand()- 0.5])
hold on
l = line([0.5 3.5 ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
title('Hyp 1 - Motor preparation ')
%ylabel('Classification')
ylim([ 0 100]);
xlim([ 0.5 2.5]);



subplot(2,2,2)
bar([NotChance, ChanceLevel])
hold on
l = line([0.5 2.5],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
title('Hyp 2 - Grasp Semantics ')
ylabel('Classification')
ylim([ 0 100]);
xlim([ 0.5 2.5]);

%xticklabels({'Image', 'Auditory', 'Written'})
xticklabels({ 'Auditory', 'Written'})


subplot(2,2,3)
bar([21 19; 86 74])
hold on
l = line([0.5 3.5 ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
title('Ex 3 A -  Word processing ')
ylabel('Classification')
ylim([ 0 100]);
xlim([ 0.5 2.5]);

%xticklabels({'Image', 'Auditory', 'Written'})
xticklabels({ 'Auditory', 'Written'})

subplot(2,2,4)
bar([86 74;51 47])
hold on
l = line([0.5 3.5 ],[chance_level,chance_level],'Color','r','LineStyle','--','Linewidth', 0.75);
title('Hyp 3 B - "language" processing')
%ylabel('Classification')
ylim([ 0 100]);
xlim([ 0.5 2.5]);

%xticklabels({'Image', 'Auditory', 'Written'})
xticklabels({ 'Auditory', 'Written'})



%xticklabels({'Image', 'Auditory', 'Written'})

legend('Grasps', 'PseudoWords', 'Chance Level', 'fontsize', 12)

sgtitle('SMG Cue Tuning')



