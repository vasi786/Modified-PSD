%% Importing the data

A = importdata('solution.log'); % import the log file generated while running PSD code.
x = A(:,1);
y = A(:,2);
z = A(:,3);
k = A(:,4);

%% Initialization

x1 = zeros(1,length(x));
x2 = zeros(1,length(x));
y1 = zeros(1,length(y));
y2 = zeros(1,length(y));
z1 = zeros(1,length(z));
z2 = zeros(1,length(z));

%% Looping to segregate the Pore sizes

for i = 1:length(x)
    if k(i) < 12 % give your desired value in ang
        x1(i) = x(i);
        y1(i) = y(i);
        z1(i) = z(i);
    end
    if k(i) > 12 % give your desired value in ang
        x2(i) = x(i);
        y2(i) = y(i);
        z2(i) = z(i);
    end
end

%% For red and blue representation

scatter3(x1,y1,z1,'b.');  % Blue is for Pore size of less than 4 ang [line no 9]
hold on;
scatter3(x2,y2,z2,'r.'); % Red is for pore size greater than 4 ang [line no 14]

%%  For colormap

% S = 50;
% scatter3(x,y,z,S,k,'filled');

%% Labelling the axis

xlabel('x');
ylabel('y');
zlabel('z');