function [] = random_message(p)
    do_it = randi(100);
    which_quote = randi(75966);
    if do_it > (100 - p)
        M = readtable('Quotes.csv');
        to_cell = M(which_quote, 1);
        to_cell = to_cell.(1);
        f = msgbox(to_cell{1});
    end
end