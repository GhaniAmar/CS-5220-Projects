{-# LANGUAGE DoRec #-}
{-# LANGUAGE BangPatterns #-}
module Keywords
where

import Types
import Control.Monad.Error
import Control.Exception (evaluate)
import Control.Concurrent.Spawn

-- (apply h t) evaluates t to a list t' and evaluates (h t') -------------------
applyKeyword eval env l@(List [Atom "apply", fun, args]) =
  do argl <- eval args env
     case argl of
       List r -> eval (List (fun:r)) env
       _      -> malformed l
applyKeyword _ _ l = malformed l

-- (car l) evalutes l and returns its head -------------------------------------
carKeyword eval env l@(List [Atom "car", lst]) =
  do list <- eval lst env
     case list of
       List (h:_)     -> return h
       Dotted (h:_) _ -> return h
       List []        -> throwError (EmptyList l)
       _              -> throwError (TypeError l)
carKeyword _ _ l = malformed l

-- (cdr l) evalutes l and returns its tail -------------------------------------
cdrKeyword eval env l@(List [Atom "cdr", lst]) =
  do list <- eval lst env
     case list of
       List (_:t)     -> return (List t)
       Dotted [] d    -> return d
       Dotted (_:t) d -> return $ Dotted t d
       List []        -> throwError (EmptyList l)
       _              -> throwError (TypeError l)
cdrKeyword _ _ l = malformed l

-- (cond pair1 pair2 ...) behaves in the usual Scheme fasion -------------------
condKeyword eval env l@(List (Atom "cond" : pairs)) = step pairs
  where step ((List [Atom "else", expr]) : _) = eval expr env
        step [] = return $ List []
        step ((List [pred, conseq]) : t) =
          do pred' <- toBool =<< eval pred env
             if pred'
               then eval conseq env
               else step t
        step _ = malformed l
condKeyword _ _ l = malformed l

-- (cons h t) evalutes h and t and forms the cons pair (h t) -------------------
consKeyword eval env (List [Atom "cons", h, t]) =
  do head <- eval h env
     tail <- eval t env
     return $ case tail of
       List l     -> List $ head : l
       Dotted l d -> Dotted (head : l) d
       _          -> Dotted [head] tail
consKeyword _ _ l = malformed l

-- (display expr) evalutes expr and prints its value to the terminal -----------
displayKeyword eval env (List [Atom "display", val]) =
  do expr <- eval val env
     liftIO (print expr)
     return (List [])
displayKeyword _ _ l = malformed l

-- (eq? a b) returns true if a and b are identical primitives or symbols -------
eq e1 e2 = case (e1, e2) of
  (Atom s1, Atom s2)           -> s1 == s2
  (Primitive p1, Primitive p2) -> p1 == p2
  _                            -> False

eqKeyword eval env (List [Atom "eq?", e1, e2]) =
  do e1' <- eval e1 env
     e2' <- eval e2 env
     return . Primitive . Bool $ eq e1' e2'
eqKeyword _ _ l = malformed l

-- (equal? a b) is like eq? but with pairwise semantics on lists ---------------
equalKeyword eval env (List [Atom "equal?", e1, e2]) =
  do e1' <- eval e1 env
     e2' <- eval e2 env
     return . Primitive . Bool $ case (e1', e2') of
       (Atom s1, Atom s2)           -> s1 == s2
       (Primitive p1, Primitive p2) -> p1 == p2
       (List l1, List l2)           -> step l1 l2
       (Dotted d1 e1, Dotted d2 e2) -> eq e1 e2 && step d1 d2
       (_, _)                       -> False
  where step [] []           = True
        step (h1:t1) (h2:t2) = eq h1 h2 && step t1 t2
        step _ _             = False
equalKeyword _ _ l = malformed l

-- (eval expr) evaluates expr twice and returns its value ----------------------
evalKeyword eval env (List [Atom "eval", expr]) =
  do expr' <- eval expr env
     eval expr' env
evalKeyword _ _ l = malformed l

-- (if pred conseq alt) is the usual if ----------------------------------------
ifKeyword eval env (List [Atom "if", pred, conseq, alt]) =
  do pred' <- toBool =<< eval pred env
     if pred'
        then eval conseq env
        else eval alt env
ifKeyword _ _ l = malformed l

-- (lambda (arg1 arg2 ...) body) does not allow zero argument functions --------
lambdaKeyword _ env (List [Atom "lambda", (List ((Atom arg) : args)), body]) =
  return $ case args of
    [] -> Closure body arg env
    _  -> Closure (List [Atom "lambda", List args, body]) arg env
lambdaKeyword _ _ l = malformed l

-- The let keyword is sequential and provides sugar for function bindings ------
process (List l) = mapM desugar l
  where desugar (List [Atom id, binding]) = return (id, binding)
        desugar (List [List ((Atom id) : args), binding]) =
          return (id, List [Atom "lambda", List args, binding])
        desugar _ = throwError . Malformed . List $ l
process _ = undefined

letKeyword eval env (List [Atom "let", bindings, body]) =
  do pairs <- process bindings
     env' <- foldM evalOne env pairs
     eval body env'
  where evalOne env (id, binding) =
          do result <- eval binding env
             return $ bind id result env
letKeyword _ _ l = malformed l

-- Only recursive functions are allowed, but mutual recursion is permitted -----
letrecKeyword eval env (List [Atom "letrec", bindings, body]) =
  do pairs <- process bindings
     rec tied <-  mapM (tiePair newEnv) pairs
         newEnv <- return $ foldl insert env tied
     eval body newEnv
  where tiePair env (id, binding) =
          do result <- eval binding env
             return (id, result)
        insert env (id, binding) = bind id binding env
letrecKeyword _ _ l = malformed l

-- (list a b ...) evaluates its arguments and returns them as a list -----------
listKeyword eval env (List (Atom "list" : tail)) =
  do elems <- mapM (`eval` env) tail
     return $ List elems
listKeyword _ _ l = malformed l

-- (list? expr) evalutes expr and returns true if it is list/dotted ------------
listqKeyword eval env (List [Atom "list?", a]) =
  do a' <- eval a env
     return . Primitive . Bool $ case a of
       List _     -> True
       Dotted _ _ -> True
       _          -> False
listqKeyword _ _ l = malformed l

-- (null? expr) returns true if and only if expr evalutes to an empty list -----
nullKeyword eval env (List [Atom "null?", a]) =
  do a' <- eval a env
     return . Primitive . Bool $ case a' of
       List [] ->  True
       _       ->  False
nullKeyword _ _ l = malformed l

-- Experimental parallel let support  ------------------------------------------
ioEval eval env expr = do expr' <- runLisp $ eval expr env
                          val <- evaluate expr'
                          case val of
                            Left err  -> error "Parallel evaluation failure"
                            Right val -> return val

manyIOEval eval exprs env = parMapIO (ioEval eval env) exprs
manyWrapper eval exprs env =
  do vals <- liftIO $ manyIOEval eval exprs env
     return vals

pletKeyword eval env l@(List [Atom "parallel-let", bindings, body]) =
  do pairs <- process bindings
     !pairs' <- manyWrapper eval (map snd pairs) env
     let env' = foldl (\env (id, expr) -> bind id expr env) env (zip (map fst pairs) pairs')
     eval body env'
  where evalOne (id, expr) =
          do expr' <- eval expr env
             return (id, expr')

-- Experimental parallel map support -------------------------------------------
pMapKeyword eval env o@(List [Atom "parallel-map", f, l]) =
  do list <- eval l env
     function <- eval f env
     case list of
       List lst -> do let jobs = map (\x -> List [function, x]) lst
                      !results <- manyWrapper eval jobs env
                      return . List $ results
       _ -> throwError (TypeError o)

-- (quote expr) returns expr untouched -----------------------------------------
quoteKeyword _ _ (List [Atom "quote", expr]) = return expr
quoteKeyword _ _ l = malformed l