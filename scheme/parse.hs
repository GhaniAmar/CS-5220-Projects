module Parse (parseString, parseFile, parseLine)
where

import Types
import Text.ParserCombinators.Parsec
import Text.Parsec.Token
import Text.Parsec.Language
import Control.Monad.Error
import System.IO

start :: Parser Char
start = oneOf "!$%&|*+-/:<=>?@^~#"

-- Lexers
scheme :: LanguageDef ()
scheme = emptyDef
    {commentStart = "",
     commentEnd = "",
     commentLine = ";",
     nestedComments = False,
     identStart = letter <|> start,
     identLetter = letter <|> start <|> char '_',
     caseSensitive = False}

lexer = makeTokenParser scheme

parsePrim :: Parser LispExpr
parsePrim =
      plift String stringLiteral
  <|> try (do
        char '-'
        plift (Double . (\x -> -x)) float)
  <|> try (plift Double float)
  <|> plift Integer integer
    where plift c l = liftM (Primitive . c) $ l lexer

parseExpr :: Parser LispExpr
parseExpr =
      liftM Atom (identifier lexer)
  <|> parsePrim
  <|> do
        char '\''
        x <- parseExpr
        return $ List [Atom "quote", x]
  <|> try list
  <|> dotted
  where list = parens lexer (liftM List $ sepBy parseExpr spaces)
        dotted = parens lexer (do head <- endBy parseExpr spaces
                                  tail <- char '.' >> spaces >> parseExpr
                                  return $ Dotted head tail)

parseString :: String -> Either ParseError [LispExpr]
parseString = parse (many (whiteSpace lexer >> parseExpr)) "scheme522"

parseLine :: String -> LispMonad [LispExpr]
parseLine l = case parseString l of
  Left err -> throwError (SyntaxError err)
  Right v  -> return v

parseFile :: FilePath ->  LispMonad [LispExpr]
parseFile p = do
  contents <- liftIO $ readFile p
  case parseString contents of
    Left err -> throwError (SyntaxError err)
    Right v  -> return v