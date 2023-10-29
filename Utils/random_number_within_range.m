function [val] = random_number_within_range(min, max, count)
    val = min+rand(1,count).*(max-min);
end