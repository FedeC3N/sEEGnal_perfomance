function selected_user = validate_user()
    
    % Read the users
    fid = fopen('users.txt','rt');
    tmp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    users = tmp{1};

    fprintf('Elige un usuario \n')
    fprintf('%s - ',users{:}) 
    fprintf('\n')
    
    % Ask for the user
    selected_user = input('Quién eres: ', 's');
    
    while true
        
        % Check the selection
        if ismember(selected_user,users)
            break
        else
            clc
            fprintf('%s no está en la lista.\n',selected_user)
            fprintf('%s - ',users{:})
            fprintf('\n')
            selected_user = input('Elegir otra vez: ', 's');
            fprintf('\n')
            
        end
        
    end
    
    fprintf('\n')
    clc

    
end