# Solitaire Texas Hold'em Poker

A single-player, graphical take on Texas Hold'em, built in Python with an object-oriented design. You play against an automated dealer, deciding after each round of cards whether to **stay** in the hand or **fold** — betting not with chips, but with points that reward good judgment about when to hold and when to walk away.

Built by **Subodh Niroula** as a computer science class project, with a focus on applying Object-Oriented Programming (OOP) principles to a real, interactive application.

## How to run it

**Requirements:** Python 3, and the `graphics.py` module included in this repo (it wraps Tkinter, which ships with most standard Python installations).

1. Make sure `graphics.py` is in the same folder as the main game script.
2. From that folder, run:
   ```bash
   python3 poker.py
   ```
   (replace `poker.py` with the actual filename of the main script if you've named it differently)
3. A game window will open. Play using your mouse — there's no keyboard input needed.

> **Note for Linux users:** if you get an error about Tkinter not being found, install it with `sudo apt-get install python3-tk` (Debian/Ubuntu) or your distribution's equivalent, then try again.

## How the game works

Each round deals out community cards in stages, and after each stage you choose to **stay** or **fold**:

1. **Hole cards.** You and the dealer each get two cards. Yours are face up; the dealer's are face down. You decide: stay or fold?
2. **The flop.** Three more community cards are dealt face up. Stay or fold?
3. **The turn.** One more community card is dealt. Stay or fold?
4. **The river.** The final community card is dealt. Stay or fold?

If you stay through all four rounds, the dealer's hole cards are revealed and the two five-card hands are compared using standard poker hand rankings. Best hand wins; a tie means no one scores.

### Scoring

**If you stay the whole way:**
- Win → **+100 points**
- Lose → **−100 points**

**If you fold early:** the game reveals what would have happened anyway, so you learn whether folding was the right call.
- If you *would have lost* — folding pays off, and the earlier you folded, the more it's worth:
  | Folded on | Points |
  |---|---|
  | Round 1 (hole cards) | +100 |
  | Round 2 (flop) | +75 |
  | Round 3 (turn) | +50 |
  | Round 4 (river) | +25 |
- If you *would have won* — folding costs you, and the later you folded, the smaller the loss:
  | Folded on | Points |
  |---|---|
  | Round 4 (river) | −25 |
  | Round 3 (turn) | −50 |
  | Round 2 (flop) | −75 |
  | Round 1 (hole cards) | −100 |

Your running **average score** is displayed after every hand, so you can track how well you're reading the game over time. After each hand, the deck is reshuffled and you're offered a new one — or you can quit at any point.

## Project structure

| File | Purpose |
|---|---|
| `poker.py` *(or your main script's actual name)* | Game logic: dealing, hand evaluation, scoring, and the UI flow |
| `graphics.py` | A simple object-oriented graphics library (by John Zelle, bundled with the *Python Programming: An Introduction to Computer Science* textbook) that wraps Tkinter for drawing windows, shapes, text, and buttons |

Under the hood, the game logic is split across several classes: `Card` and `Hands` for dealing, `CountingHands` and `ComparingHands` for evaluating poker hands (straight flush down to high card), `DrawCards` and `Button` for rendering the table, and `StartGame` for orchestrating a full round from deal to payout.

## Known limitations

- The window is a fixed size and isn't resizable.
- All interaction is mouse-based (click the on-screen buttons); there's no keyboard shortcut support.
- Hand comparison in a rare case (both players landing exactly on "One Pair") falls back to comparing the pair's rank directly rather than the full standard poker tiebreaker rules (e.g., kickers aren't considered).

## Acknowledgment

Thank you to Professor Anna Varvak for assigning this project and for the guidance and support throughout its development.
