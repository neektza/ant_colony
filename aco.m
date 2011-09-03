clearvars;

%% init values
load('graph.mat', 'graph');

ro = 0.2; % pridonos najboljeg puta
fi = 0.1; % zaboravljanje

len = length(graph); % duzina grafa

tao = zeros(size(graph)); % feromoni grane (trail level)
ni = 1./graph; % privlacnost grane (attractiveness)

s = len/2; % broj grana grafa u kojima odlucujemo nasumicno
q = 1 - s/len; % probability of exploitation

t_max = 500; % max broj iteracija
interval_max = t_max/4; % bax broj iteracija bez poboljsanja
sfl = max(graph(:)); % shortest found length
ants = len; % broj mravaca

% inicjalizacija matrice vjerojatnosti
probs = zeros(1,len);

%% inital trail level calculation
paths = zeros(ants,len); % init/reset ruta mravaca
paths(:,1) = 1; % svi mravci krecu iz cvora 1
paths(:,27) = 27; % svi mravci zavrsavaju u cvoru 27

for ant = 1 : ants
    g = graph; % radna kopija grafa (za svakog mrava posebna)
    for s = 1 : len-1 % u 26. cvoru je zadnja odluka

    cn = paths(ant, s); % cn, current node
    g(g(:,cn) == -1 ,cn) = -2; % u kolumni posjecenog cvora, oznaci da vise nije zabranjen za taj cvor
    
    vn = paths(ant, 1:s); % vn, visited nodes   
    [r c]  = find(g == -1); % forbidden nodes - svi cvorovi koji u redu imaju barem jednu jedinicu
    fn = transpose(unique(r));

    probs = ni(cn, 1:len);
    probs(fn) = 0; % ne mozemo u -1 cvorove
    probs(vn) = 0; % ne mozemo u one u kojima smo bili

    eve = randsample([1, 0], 1, true, [q, 1-q]);

    if eve % best
        i = find(probs == max(probs));
        nn = i(randi(length(i),1,1));
    else % random
        nn = randsample(1:len, 1, true, probs);
    end

    paths(ant, s+1) = nn; % nn, next node
    end
end

% izracunaj duzine tura
path_lengths = zeros(size(paths,1),1);
for ant = 1 : size(paths,1) 
    for s = 1 : len-1
        path_lengths(ant) = path_lengths(ant) + graph(paths(ant,s),paths(ant,s+1));
    end
end

% postavi inicjalni trail level, tao = (Lnn * n)^-1
best_length = min(path_lengths);
t0 = (best_length * len) ^ -1;

%% acs - sop

t = 1; timestamp = 1; tao(:) = t0;
while((t <= t_max) && (t-timestamp < interval_max))
    
    paths = zeros(ants,len); % init/reset ruta mravaca
    paths(:,1) = 1; % svi mravci krecu iz cvora 1
    paths(:,27) = 27; % svi mravci zavrsavaju u cvoru 27

    for ant = 1 : ants
        g = graph; % radna kopija grafa (za svakog mrava posebna)
        for s = 1 : len-1 % u 26. cvoru je zadnja odluka

        cn = paths(ant, s); % cn - current node
        g(g(:,cn) == -1 ,cn) = -2; % u kolumni posjecenog cvora, oznaci da vise nije zabranjen za taj cvor

        vn = paths(ant, 1:s); % vn - visited nodes   
        [r c]  = find(g == -1); % forbidden nodes - svi cvorovi koji u redu imaju barem jednu jedinicu
        fn = transpose(unique(r));

        probs = tao(cn,1:len) .* ni(cn, 1:len);
        probs(fn) = 0; % ne mozemo u -1 cvorove
        probs(vn) = 0; % ne mozemo u one u kojima smo bili

        eve = randsample([1, 0], 1, true, [q, 1-q]); % e vs e

        if eve % best
            i = find(probs == max(probs));
            nn = i(randi(length(i),1,1));
        else % random
            nn = randsample(1:len, 1, true, probs);
        end
        
        tao(cn,nn) = (1-fi) * tao(cn,nn) + fi * t0;

        paths(ant, s+1) = nn; % nn - next node
        end
    end

    % duzine puteva
    path_lengths = zeros(size(paths,1),1);
    for ant = 1 : size(paths,1) 
        for s = 1 : len-1
            path_lengths(ant) = path_lengths(ant) + graph(paths(ant,s),paths(ant,s+1));
        end
    end

    best_length = min(path_lengths);
    best_ant = find(path_lengths == best_length);
    best_path = paths(best_ant, :);

    for s = 1 : len-1 % for all nodes along the best path, place feromones on edges
        i = best_path(s);
        j = best_path(s+1);
        tao(i,j) = (1-ro)*tao(i,j) + ro/best_length;
    end
    
    if best_length < sfl
        sfl = best_length;
        sfp = best_path;
        timestamp = t; % iteracija zadnjeg poboljsanja
    end
    
    t = t+1;
    best_length
end

sfl
sfp

