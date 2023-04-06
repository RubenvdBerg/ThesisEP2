from itertools import combinations

game_info = {'Ruben': 'Game A', 'Marlo': 'Game B', 'Elmer': 'Game C', 'Eric': 'Game D', 'Daan': 'Game E', 'Joost': 'Game F', 'Jan': 'Game G'}
names = game_info.keys()
games = game_info.values()

match_ups = list(combinations(names, 2))

# Dictionary to keep track of which game each individual has been assigned
games_played = {name: set() for name in names}

game_list = dict.fromkeys(match_ups)
# Choose the game each name combination will play
for name_combo in match_ups:
    for game in games:
        # Check if either of the two people has not already played this game
        if not any(game in games_played[name] for name in name_combo):
            # Add the game to games played of both people
            for name in name_combo:
                games_played[name].add(game)
            # Assign the game to this name combination
            game_list[name_combo] = game
            break
# Check if every matchup got assigned a game
if None in game_list.values():
    print(game_list)
    raise ValueError('Not all match-ups have been assigned a game')

rounds = {i: [] for i in range(1, len(game_info) + 1)}
round_sets = {i: set() for i in rounds}
for (player1, player2), game in game_list.items():
    contest = [player1, player2, game]
    for round_nr in rounds:
        if contest not in rounds[round_nr]:
            if not any(x in round_sets[round_nr] for x in contest):
                rounds[round_nr].append(contest)
                round_sets[round_nr].update(contest)
                break

for round_nr, games in rounds.items():
    print(f'Round #{round_nr}')
    for i, (player1, player2, game_type) in enumerate(games):
        print(f'Game {i + 1}: {player1:>6}, {player2:>6}, {game_type}')
