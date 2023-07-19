function nameOut = simplifyName(nameIn)
% written Jan 2023 by M. Osman, Cambridge University (mo549@cam.ac.uk)
% tokenizes an input name to a lower case first initial and lower case surname

if ~ismissing(nameIn)

nameIn = lower(nameIn);
nameIn = regexprep(nameIn,'(?:á|à|â|ä|ã|å)','a');
nameIn = regexprep(nameIn,'(?:æ)','ae');
nameIn = regexprep(nameIn,'(?:ç)','c');
nameIn = regexprep(nameIn,'(?:ð)','d');
nameIn = regexprep(nameIn,'(?:é|è|ê|ë)','e');
nameIn = regexprep(nameIn,'(?:í|ì|î|ï)','i');
nameIn = regexprep(nameIn,'(?:ñ)','n');
nameIn = regexprep(nameIn,'(?:ó|ò|ô|ö|õ|ø)','o');
nameIn = regexprep(nameIn,'(?:œ)','oe');
nameIn = regexprep(nameIn,'(?:ú|ù|ü|û)','u');
nameIn = regexprep(nameIn,'(?:ý|ÿ)','y');
nameIn = regexprep(nameIn,'(?:š)','s');
nameIn = regexprep(nameIn,'(?:č)','c');
% replace special characters
nameIn = replace(nameIn,', ','');
nameIn = replace(nameIn,',',''); 
nameIn = replace(nameIn,'.','');
nameIn = replace(nameIn,"'",'');

cp = nameIn; 
cp = split(cp); 
if length(char(cp(1))) > 1
    firstletter = char(cp(1)); 
    cp(1) = string(firstletter(1)); 
end
nameOut = strcat(cp(1)," ",cp(end)); 

else

    warning('Input name is missing!')
    nameOut = nameIn; 

end

end