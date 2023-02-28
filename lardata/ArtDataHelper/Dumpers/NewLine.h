/**
 * @file   NewLine.h
 * @brief  Simple class managing a repetitive output task
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 18th, 2015
 */

#ifndef LARDATA_RECOBASE_DUMPERS_NEWLINE_H
#define LARDATA_RECOBASE_DUMPERS_NEWLINE_H 1

// C/C++ standard libraries
#include <string>
#include <utility> // std::move()

namespace recob {
  namespace dumper {

    /// Structure collecting indentation options
    struct IndentOptions_t {
      std::string indent;       ///< indentation string
      bool appendFirst = false; ///< skip indentation on the first line

      IndentOptions_t(std::string ind = "", bool followLine = false)
        : indent(ind), appendFirst(followLine)
      {}

      IndentOptions_t& appendIndentation(std::string more)
      {
        indent += more;
        appendFirst = false;
        return *this;
      }
      IndentOptions_t& removeIndentation(std::string less)
      {
        indent.erase(std::max(indent.length() - less.length(), size_t(0)));
        return *this;
      }

    }; // IndentOptions_t

    /**
     * @brief Starts a new line in a output stream
     * @tparam Stream type of output stream
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * std::cout << "Preamble on its own line." << std::endl;
     * NewLine OutLn(std::cout, "> ");
     * OutLn() << "An indented line.";
     * OutLn() << "Another indented line.";
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * that (after flush) will result in the output
     *
     *     Preamble on its own line.
     *
     *     > An indented line.
     *     > Another indented line.
     *
     * Asking to consider the first line already started, instead:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * std::cout << "Preamble on its own line." << std::endl;
     * NewLine OutLn(std::cout, "> ", true);
     * OutLn() << "An indented line.";
     * OutLn() << "Another indented line.";
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will instead result in the output
     *
     *     Preamble on its own line.
     *     > An indented line.
     *     > Another indented line.
     *
     * (note that the line that we consider started was actually an empty one).
     */
    template <typename Stream>
    class NewLine {
    public:
      /**
       * @brief Constructor: associates with the stream
       * @param stream a reference to the stream where to insert new lines
       * @param indentOptions all indentation options (will be copied)
       *
       * The constructor does not start a new line.
       * If followLine is true, the first line is supposed to be already started
       * and no indentation nor new line will be set on it.
       */
      NewLine(Stream& stream, IndentOptions_t indentOptions)
        : out(stream), options(std::move(indentOptions)), nLines(0)
      {}

      /**
       * @brief Constructor: associates with the stream
       * @param stream a reference to the stream where to insert new lines
       * @param indent string used for indentation (default: none)
       * @param followLine whether first line is already started (default: no)
       *
       * The constructor does not start a new line.
       * If followLine is true, the first line is supposed to be already started
       * and no indentation nor new line will be set on it.
       */
      NewLine(Stream& stream, std::string indent = "", bool followLine = false)
        : NewLine(stream, IndentOptions_t{indent, followLine})
      {}

      /// @{
      /// @name Accessors

      /// Returns the number of inserted lines
      unsigned int lines() const { return nLines; }

      /// Returns the current indentation string
      std::string indent() const { return options.indent; }

      /// @}

      /// Starts a new line
      Stream& newLine()
      {
        if (!append()) forceNewLine();
        ++nLines;
        return out;
      }

      /// Calls and returns newLine(). Candy.
      Stream& operator()() { return newLine(); }

      /// Starts a new line (no matter what)
      void forceNewLine() { out << "\n" << options.indent; }

      /// Returns whether newLine() will append text on the current line
      bool append() const { return (lines() == 0) && options.appendFirst; }

      /// Replaces the indentation string
      void setIndent(std::string newIndent) { options.indent = newIndent; }

      /// Adds to the end to the indentation string
      void addIndent(std::string moreIndent) { options.indent += moreIndent; }

    protected:
      Stream& out;             ///< reference to the output stream
      IndentOptions_t options; ///< all indentation options
      unsigned int nLines;     ///< number of lines in output

    }; // class NewLine

    /// Convenience function to create a temporary NewLine
    template <typename Stream>
    inline NewLine<Stream> makeNewLine(Stream& stream, std::string indent, bool followLine = false)
    {
      return NewLine<Stream>(stream, indent, followLine);
    }

    /// Convenience function to create a temporary NewLine
    template <typename Stream>
    inline NewLine<Stream> makeNewLine(Stream& stream, IndentOptions_t const& options)
    {
      return NewLine<Stream>(stream, options);
    }

  } // namespace dumper
} // namespace recob

#endif // LARDATA_RECOBASE_DUMPERS_NEWLINE_H
