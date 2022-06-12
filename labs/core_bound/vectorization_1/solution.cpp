#include "solution.hpp"
#include <cassert>
#include <type_traits>

using vectorized_val = std::array<uint8_t, sequence_count_v>;
using vectorized_sequence_T = std::array<vectorized_val, sequence_size_v>;

vectorized_sequence_T convert_to_seq_matrix_T (std::vector<sequence_t> const &sequences1) {
  vectorized_sequence_T seq_mat;
  for (int row = 0; row < sequence_size_v; ++row) {
    for (int col = 0; col < sequence_count_v; ++col) {
      seq_mat[row][col] = sequences1[col][row];
    }
  }

  return seq_mat;
}

result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
    auto vectorized_seq_1 = convert_to_seq_matrix_T(sequences1);
    auto vectorized_seq_2 = convert_to_seq_matrix_T(sequences2);

    using score_t = int16_t;
    using vectorized_score_t = std::array<score_t, sequence_count_v>;
    using vectorized_score_column_t = std::array<vectorized_score_t, sequence_size_v + 1>;
    /*
     * Initialise score values.
     */
    score_t gap_open{-11};
    score_t gap_extension{-1};
    score_t match{6};
    score_t mismatch{-4};

    /*
     * Setup the matrix.
     * Note we can compute the entire matrix with just one column in memory,
     * since we are only interested in the last value of the last column in the
     * score matrix.
     */
    vectorized_score_column_t score_column{};
    vectorized_score_column_t horizontal_gap_column{};
    vectorized_score_t last_vertical_gaps{};

    /*
     * Initialise the first column of the matrix.
     */
    horizontal_gap_column[0].fill(gap_open);
    last_vertical_gaps.fill(gap_open);

    for (size_t i = 1; i < score_column.size(); ++i) {
      for (size_t j = 0; j < sequence_count_v; ++j) {
        score_column[i][j] = last_vertical_gaps[j];
        horizontal_gap_column[i][j] = last_vertical_gaps[j] + gap_open;
        last_vertical_gaps[j] += gap_extension;
      }
    }

    /*
     * Compute the main recursion to fill the matrix.
     */
    for (unsigned col = 1; col <= vectorized_seq_2.size(); ++col) {
      auto last_diagonal_score = score_column[0]; // Cache last diagonal score to compute this cell.
      score_column[0] = horizontal_gap_column[0];
      for (int i = 0; i < sequence_count_v; ++i) {
        last_vertical_gaps[i] = horizontal_gap_column[0][i] + gap_open;
        horizontal_gap_column[0][i] += gap_extension;
      }      

      for (unsigned row = 1; row <= vectorized_seq_1.size(); ++row) {
        // Compute next score from diagonal direction with match/mismatch.
        vectorized_score_t best_cell_score = last_diagonal_score;
        for (int i = 0; i < sequence_count_v; ++i) {
          best_cell_score[i] +=
              (vectorized_seq_1[row - 1][i] == vectorized_seq_2[col - 1][i] ? match : mismatch);
        }

        for (int i = 0; i < sequence_count_v; ++i) {
          // Determine best score from diagonal, vertical, or horizontal
          // direction.
          best_cell_score[i] = max(best_cell_score[i], last_vertical_gaps[i]);
          best_cell_score[i] = max(best_cell_score[i], horizontal_gap_column[row][i]);
          // Cache next diagonal value and store optimum in score_column.
          last_diagonal_score[i] = score_column[row][i];
          score_column[row][i] = best_cell_score[i];
          // Compute the next values for vertical and horizontal gap.
          best_cell_score[i] += gap_open;
          last_vertical_gaps[i] += gap_extension;
          horizontal_gap_column[row][i] += gap_extension;
          // Store optimum between gap open and gap extension.
          last_vertical_gaps[i] = max(last_vertical_gaps[i], best_cell_score[i]);
          horizontal_gap_column[row][i] =
              max(horizontal_gap_column[row][i], best_cell_score[i]);
        }
      }
    }

    return score_column.back();
}