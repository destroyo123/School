function ua = uacolor()
%% Primary
ua.red = uint8([171,5,32]);
ua.blue = uint8([12,35,75]);

%% Primary complimentary
% Reds
ua.bloom = uint8([239 64 86]);
ua.chili = uint8([139 0 21]);

% Blue
ua.sky = uint8([129 211 235]);
ua.oasis = uint8([55 141 189]);
ua.azurite = uint8([30 82 136]);
ua.midnight = uint8([0 28 72]);

%% Neutral colors
ua.cool_gray = uint8([226 233 235]);
ua.warm_gray = uint8([244 237 229]);

%% Secondary colors
ua.leaf = uint8([112 184 101]);
ua.river = uint8([0 125 132]);
ua.silver = uint8([158 171 174]);
ua.mesa = uint8([169 92 66]);

%% Legacy colors
ua.ash = uint8([64 54 53]);
ua.sage = uint8([74 99 78]);
ua.copper = ua.mesa;

end