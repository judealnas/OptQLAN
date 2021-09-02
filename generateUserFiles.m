user_arr = [];
for i = 1:N
    user_arr = [user_arr ; struct('name', user_names(i),...
                                  'eff', eff_N(i),...
                                  'dark_count_rate', d_N(i))];
end

disp(user_arr)
json = jsonencode(user_arr, 'PrettyPrint', true);
disp(json)

user_file = fopen('users.txt','w');
fprintf(user_file, "%s", json);
fclose(user_file);