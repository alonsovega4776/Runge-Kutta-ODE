%%Data analysis 
close all
clear

%% Preproccesing 
file =  fopen('./RungeKuttaODE_C11/cmake-build-debug/output.txt', 'r');
if (file == 0)
    fprintf('ERROR: Can not open file.\n')
    exit
end
fprintf('Opened file. \n')


t_str = fgetl(file);
t_str = t_str(find(t_str == 32, 1) + 1:find(t_str == 44, 1, 'last') - 1);
t_str = convertCharsToStrings(t_str);
t = cell2mat(textscan(t_str,'%f', 'Delimiter',','));

x_1_str = fgetl(file);
x_1_str = x_1_str(find(x_1_str == 32, 1) + 1:find(x_1_str == 44, 1, 'last') - 1);
x_1_str = convertCharsToStrings(x_1_str);
x_1 = cell2mat(textscan(x_1_str,'%f', 'Delimiter',','));

x_2_str = fgetl(file);
x_2_str = x_2_str(find(x_2_str == 32, 1) + 1:find(x_2_str == 44, 1, 'last') - 1);
x_2_str = convertCharsToStrings(x_2_str);
x_2 = cell2mat(textscan(x_2_str,'%f', 'Delimiter',','));

%% Plotting 
figure
plot(t,x_1);
hold on
hold off

figure
plot(t,x_2);
hold on
hold off

%% Close files

fclose(file);